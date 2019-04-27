#===============================================================================
#                           BETA2 DMS PIPELINE
#                           Nathan Lubock
#===============================================================================

SHELL := /bin/bash
python2 := /opt/conda/envs/py27/bin/python
python3 := python

MAPTHREADS ?= 4 # how many threads to run the barcode mapper (~60 Gb per thread for all reads)
THREADS ?= 10 # how many threads to run the BB* portion of pipeline
COMPTHREADS ?= 4 # how many compression threads

all: map qseqs fastqs neg-controls synon cons

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# map-based recipes

MAPRUNS := $(addprefix pipeline/, $(addsuffix .merge.fastq, \
    Eric_S5 Map_Lane_1 Map_Lane_2 Map_Lane_3 Map_Lane_4 \
    Map-2_Lane_1 Map-2_Lane_2 Map-2_Lane_3 Map-2_Lane_4))

map: output/NextSeq_MiSeq.known-vars.txt.gz
neg-controls: output/NextSeq_MiSeq.negs.txt.gz
synon: output/NextSeq_MiSeq.synon.translate.txt

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# rna-seq recipes

QSEQS := drug-2 drug-3
FASTQS := drug-6 drug-7
qseqs: $(addprefix output/, $(addsuffix _idx-bcs-counts.txt.gz, $(QSEQS)))
fastqs: $(addprefix output/, $(addsuffix _cond-bcs-counts.txt.gz, $(FASTQS)))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# other

cons: $(addprefix ancillary/cons/, $(addsuffix _js.tsv, species class-a))

clean:
	rm -f pipeline/*

.PRECIOUS: $(addprefix pipeline/, %.map.csv %.merge.fastq %.filter.fastq \
    %_idx-bcs.txt.gz %_neg-control_vars.txt)

#===============================================================================
#                     BARCODE MAPPING
#===============================================================================

pipeline/phiX.fasta:
	@echo "Grabbing the PhiX genome"
	@curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&amp;id=NC_001422.1&amp;rettype=fasta&amp;retmode=text" >> $@

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Read Processing Pipeline:
# -------------------------
# 1) Filter out PhiX Reads
# 2) Trim adapter sequences (p5/p7/primers)
# 3) Merge reads
# 4) Filter out any reads with an N basecall
# 5) Trim any remaining reads with left-over adapters (especially 5' end)
pipeline/%.merge.fastq: pipeline/phiX.fasta data/map/%_R1.fastq.gz data/map/%_R2.fastq.gz data/seq-adapters.fasta
	@echo Filtering, merging, and trimming - $(word 2, $^)
	@bbduk2.sh \
	    in1=$(word 2, $^) \
	    in2=$(word 3, $^) \
	    fref=$< \
	    rref=$(word 4, $^) \
	    k=23 \
	    mink=11 \
	    hdist=1 \
	    trimbyoverlap=t \
	    trimpairsevenly=t \
	    overwrite=t \
	    threads=$(THREADS) \
	    stats=$(@:.merge.fastq=.filter.stats.txt) \
	    out=stdout.fastq \
	    -Xmx8g 2> $(@:.merge.fastq=.filter.err) | \
	    bbmerge.sh \
	    in=stdin.fastq \
	    interleaved=t \
	    threads=$(THREADS) \
	    adapters=$(word 4, $^) \
	    outm=stdout.fastq \
	    2> $(@:.fastq=.err) | \
	    bbduk2.sh \
	    in=stdin.fastq \
	    rliteral='GGTCGCCCTTATTACTACCAAGCTCGTGGACGGAGGC' \
	    lliteral='AAGTGCCTTCCTGCCCTTTAATCAGATGCGTCG' \
	    k=18 \
	    mink=11 \
	    hdist=1 \
	    maxns=0 \
	    overwrite=t \
	    threads=$(THREADS) \
	    stats=$(@:.merge.fastq=.trim.stats.txt) \
	    out=$@ \
	    -Xmx8g 2> $(@:.merge.fastq=.trim.err)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# extract and uniq rna-seq bcs
pipeline/rna-bcs.txt: $(addprefix output/, $(addsuffix -bcs-counts.txt.gz, drug-3_idx drug-6_cond drug-7_cond))
	zcat $^ | \
	    awk '{a[$$(NF - 1)]++} END {for(bc in a) print bc}' > $@

# Barcode Mapper:
# ---------------
#  1) Only consider barcodes from the RNA-seq data
#  2) Enforce a number of quality controls (e.g min number of reads, check that
#     variants map to the same ADRB2 chunk, no truncations, no variants from
#     the same chunk; see methods/bc-map.py source for more details)
#  3) Collapse on majority base-call
pipeline/NextSeq_MiSeq.map.csv: pipeline/rna-bcs.txt data/ADRB2.fasta $(MAPRUNS)
	@echo "Barcode mapping - $(filter-out $(wordlist 1, 2, $^), $^)"
	@cat $(filter-out $(wordlist 1, 2, $^), $^) | \
	    $(python2) ./scripts/bc-map.py \
	    -v \
	    -j$(MAPTHREADS) \
	    --bc-start 1 \
	    --bc-length 15 \
	    --min-reads 3 \
	    --start-dist 5 \
	    --bbmap-procs $(THREADS) \
	    --trunc-len 1 \
	    --contam-reads 2 \
	    --contam-dist 4 \
	    -b $(@:.map.csv=.bad-bcs.txt) \
	    - $< $(word 2, $^) > $@.lock \
	    2> $(@:.csv=.err) && \
	    mv $@.lock $@

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Find Designed Variants in Barcode Map:
# --------------------------------------
# 1) Take Map as reference
# 2) Map designed mutants onto sequences. NOTE these designs are missing some of
#    the plasmid backbone included in the mapping step. Since these sections are
#    clonally verified prior to mapping, it is unlikely that an indel occured in
#    those regions.
# 3) BBMap's semiperfectmode=t ensures a perfect match in to the reference while
#    allowing for longer reads
output/%.known-vars.txt.gz: pipeline/%.map.csv data/ADRB2_mutants.fasta
	@echo "Mapping single point mutants for $<"
	@awk -F, '{print ">"$$1"_"$$3"\n"$$2}' $< | \
	    bbmap.sh \
	    ref=stdin.fasta \
	    in=$(word 2, $^) \
	    outm=pipeline/$(@F:.txt.gz=.sam) \
	    nodisk=t \
	    noheader=t \
	    semiperfectmode=t \
	    maxindel=500 \
	    ambiguous=all \
	    secondary=t \
	    ssao=t \
	    maxsites=1000000 \
	    overwrite=t \
	    -Xmx64g \
	    threads=$(THREADS) 2> pipeline/$(@F:.txt.gz=.err)
	@awk '{sub(/_/, " ", $$3); print $$3, $$1}' pipeline/$(@F:.txt.gz=.sam) | \
	    pigz -c -p$(COMPTHREADS) > $@

#===============================================================================
#                        BARCODE COUNTING
#===============================================================================

# QSEQ Handling:
# --------------
# 1) ensure qseqs are named lane_cluster_read.qseq.gz
# 2) collapse into index barcode
# 3) only count barcodes from valid indices

data/qseqs/%/rename.err: data/qseqs/%/get-data.sh
	@echo "Ensuring data exists for $(<D)"
	@cd $(@D) && $(SHELL) $(<F)

pipeline/%_idx-bcs.txt.gz: data/qseqs/%/rename.err
	@echo "Grabbing barcodes for all qseqs in $(<D)"
	@parallel -j$(THREADS) --xapply bash ./scripts/qseq2txt.sh {1} {2} ::: $(<D)/*_1.qseq.gz ::: $(<D)/*_2.qseq.gz | \
	    awk '{print $$1, substr($$2, 1, 15)}' | \
	    pigz -c -p$(COMPTHREADS) > $@

output/%_idx-bcs-counts.txt.gz: pipeline/%_idx-bcs.txt.gz data/%_idx-list.txt
	@echo "Counting barcodes in $<"
	@parallel zcat $< \| \
	    awk -v name=\"{}\" \''$$1 == name {a[$$2]++} END {for(bc in a) print name, bc, a[bc]}'\' \
	    :::: $(lastword $^) | \
	    pigz -c -p$(COMPTHREADS) > $@

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# One off for drug-3:
# -------------------
# 1) collapse index read
# 2) trim off N's
# 3) filter out anything that maps to phiX
pipeline/drug-3_idx-bcs.txt.gz: data/qseqs/drug-3/rename.err pipeline/phiX.fasta
	@echo "Quality filtering reads from $(<D)"
	@parallel -j$(THREADS) --xapply bash ./scripts/qseq2fastq.sh {1} {2} ::: $(<D)/*_1.qseq.gz ::: $(<D)/*_2.qseq.gz | \
	    bbduk2.sh \
	    in=stdin.fastq \
	    outu=stdout.fasta \
	    ref=$(lastword $^) \
	    minavgquality=25 \
	    forcetrimright=14 \
	    k=14 \
	    hdist=1 \
	    overwrite=t \
	    stats=$(@:.txt.gz=.filter.txt) \
	    threads=$(THREADS) \
	    -Xmx16g  \
	    2> $(@:.txt.gz=.err) | \
	    sed 's/>\|_/ /g' | \
	    paste - - | \
	    awk '{print $$1, $$3}' | \
	    pigz -c -p$(COMPTHREADS) > $@

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# FASTQ Handling:
# ---------------
# 1) ensure files are named properly
# 2) count barcodes

data/fastqs/%/rename.err: data/fastqs/%/get-data.sh
	@echo "Ensuring data exists for $(<D)"
	@cd $(@D) && $(SHELL) $(<F)

# combine reads from 4 lanes in second awk step
output/drug-6_cond-bcs-counts.txt.gz: data/fastqs/drug-6/rename.err
	@echo "Counting barcodes for drug-6"
	@parallel -j$(THREADS) zcat {} \| \
	    awk -v name=\"{/.}\" \''NR % 4 == 2{a[substr($$1,1,15)]++} END {for(bc in a) print name, bc, a[bc]}'\' \
	    ::: $(<D)/*.fastq.gz | \
	    sed -e 's/_.\.fastq//' | \
	    awk '{a[$$1"_"$$2] += $$3} END {for(expr in a) print expr, a[expr]}' | \
	    sed -e 's/_/ /g' | \
	    pigz -c -p$(COMPTHREADS) > $@

output/drug-7_cond-bcs-counts.txt.gz: data/fastqs/drug-7/rename.err
	@echo "Counting barcodes for drug-7"
	@parallel -j$(THREADS) zcat {} \| \
	    awk -v name=\"{/.}\" \''NR % 4 == 2{a[substr($$1,1,15)]++} END {for(bc in a) print name, bc, a[bc]}'\' \
	    ::: $(<D)/*.fastq.gz | \
	    sed -e 's/.fastq//' -e 's/_/ /' | \
	    pigz -c -p$(COMPTHREADS) > $@

#===============================================================================
#                       POSITIVE/NEGATIVE CONTROLS
#===============================================================================

# Intermediate Trimming:
# ----------------------
# 1) Convert barcode map to fasta
# 2) Trim constant sequences off mapped variants
pipeline/NextSeq_MiSeq.map.trim.fasta: pipeline/NextSeq_MiSeq.map.csv data/5-prime-trim.fasta
	@echo "Trimming constant region off $<"
	@awk -F, '{print ">"$$1"_"$$3"\n"$$2}' $< | \
	    bbduk.sh \
	    in=stdin.fasta \
	    ref=$(word 2, $^) \
	    ktrim=l \
	    restrictleft=40 \
	    k=16 \
	    mink=7 \
	    edist=1 \
	    edist2=0 \
	    overwrite=t \
	    -Xmx8g \
	    threads=$(THREADS) \
	    stats=$(@:.trim.fasta=.left-trim.stats.txt) \
	    out=$@ 2> $(@:.trim.fasta=.left-trim.err)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Align map to ADRB2:
# -------------------
pipeline/NextSeq_MiSeq.map.sam: pipeline/NextSeq_MiSeq.map.trim.fasta data/ADRB2.fasta
	@echo "Aligning barcode map to ADRB2"
	@bbmap.sh \
	    in=$< \
	    ref=$(lastword $^) \
	    nodisk=t \
	    noheader=t \
	    k=8 \
	    vslow=t \
	    maxindel=500 \
	    overwrite=t \
	    threads=$(THREADS) \
	    outm=$@ \
	    outu=$(@:.sam=.no-map.sam) \
	    -Xmx16g \
	    2> $(@:.sam=.err)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Classify Negative Controls:
# ---------------------------
# 1) Grab CIGAR string from alignment
# 2) Count insertions and deletions to figure out the frame
output/NextSeq_MiSeq.negs.txt.gz: pipeline/NextSeq_MiSeq.map.sam
	@echo "Looking for frameshifts in $<"
	@awk '{print $$1, $$4, $$6}' $< | \
	    $(python2) ./scripts/classify-negs.py -j2 - | \
	    sed 's/_/\t/g' | \
	    pigz -c -p$(COMPTHREADS) > $@

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Find Perfect Variants in Barcode Map:
# -------------------------------------
pipeline/NextSeq_MiSeq.perfects.txt: pipeline/NextSeq_MiSeq.map.trim.fasta data/ADRB2.fasta
	@echo "Finding perfect alignments to ADRB2"
	@bbmap.sh \
	    in=$< \
	    ref=$(lastword $^) \
	    nodisk=t \
	    noheader=t \
	    perfectmode=t \
	    maxindel=500 \
	    overwrite=t \
	    threads=$(THREADS) \
	    outm=stdout.sam \
	    -Xmx16g 2> pipeline/$(@F:.txt=.err) | \
	    awk '{sub(/=/,""); if($$6 > 180) print $$1, $$4, $$6, $$10}' | \
	    sed 's/_/ /' > $@

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Find Synonymous Variants in Barcode Map:
# ----------------------------------------
#  1) Generate all single synonymous codon changes to the wt sequence
#  2) Align map perfectly to it
output/NextSeq_MiSeq.synon.txt: pipeline/NextSeq_MiSeq.map.trim.fasta data/ADRB2_synon.fasta
	@echo "Mapping single point mutants for $<"
	@bbmap.sh \
	    in=$< \
	    ref=$(lastword $^) \
	    nodisk=t \
	    noheader=t \
	    perfectmode=t \
	    maxindel=500 \
	    overwrite=t \
	    threads=$(THREADS) \
	    outm=stdout.sam \
	    -Xmx16g 2> pipeline/$(@F:.txt=.err) | \
	    awk '{sub(/=/,""); if($$6 > 180) print $$1, $$3, $$4, $$6, $$10}' | \
	    sed 's/_/ /' > $@

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Translate sequences to find synonymous mutations
# ------------------------------------------------
# Since our primers did not end on a codon, we need to trim certain chunks
# before translating. We use the rough bbmap alignment to get an idea of where
# each chunk aligns to. Since the leftmost part of each chunk is clonal, we
# expect this part of the alignment to be perfect. We then use synon-filter.py
# to carry out the trimming. We can use a local parasail alignment as long as we
# check that the input length equals the number of matches. We will also filter
# out anything that has a perfect nucleotide alignment to our reference as it is
# probably not actually wt
output/NextSeq_MiSeq.synon.translate.txt: pipeline/NextSeq_MiSeq.map.sam pipeline/NextSeq_MiSeq.perfects.txt data/ADRB2_prot.fasta
	@echo "Finding synonymous variants in $<"
	@$(python2) ./scripts/synon-filter.py $< | \
	    translate6frames.sh \
	    in=stdin.fasta \
	    out=pipeline/tmp.fasta \
	    tag=f \
	    overwrite=t \
	    fastawrap=10000 \
	    frames=1 \
	    2> /dev/null
	@parasail_aligner \
	    -c 60 \
	    -t $(THREADS) \
	    -f $(lastword $^) \
	    -q pipeline/tmp.fasta \
	    -a sw_stats_scan \
	    -g pipeline/$(@F:.txt=.parasail.txt)
	@# parse the index to get the barcode (parasail 0-index)
	@awk 'FNR==NR \
	    {if($$3 == $$8) a[$$1 + 1]; next} \
	    {if(FNR in a) print substr($$1,2)}' \
	    FS=',' pipeline/$(@F:.txt=.parasail.txt) \
	    FS='\t' <(paste - - < pipeline/tmp.fasta) | \
	    sed 's/_/ /g' > pipeline/$(@F:.txt=.unfilter.txt)
	@awk 'FNR==NR {a[$$1];next} !($$1 in a)' \
	    $(word 2, $^) \
	    pipeline/$(@F:.txt=.unfilter.txt) > $@

#===============================================================================
#                   ANCILLARY ANALYSIS
#===============================================================================

# convert mutational tolerance to something chimera can read
ancillary/%.4ldl.txt: ancillary/d3/%.csv
	echo "attribute: tolerance" > $@
	awk -F, 'NR > 1{print $$2, $$3}' $< | \
	    sed 's/_//' | \
	    awk '{print "\t:"$$1+1000"\t"$$2}' >> $@

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Format Jensen-Shannon scores:
# -----------------------------
# 1) ADRB2_human will be at a defined location in the JS output
# 2) Grab the conservation score associated with that position
# 3) Ignore any gaps
ancillary/cons/species_js.tsv: ancillary/cons/oma-HUMAN24043_js.raw.tsv
	awk 'NR > 2{print $$2, substr($$3, 1, 1)}' $< |\
	    awk '$$2 != "-"{print $$1, $$2}' OFS='\t' > $@

ancillary/cons/class-a_js.tsv: ancillary/cons/class-a_js.raw.tsv
	awk 'NR > 2{print $$2, substr($$3, 25, 1)}' $< |\
	    awk '$$2 != "-"{print $$1, $$2}' OFS='\t' > $@

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Grab residues between lipid membrane
ancillary/3sn6-lipid.txt: ancillary/pdb/3sn6.opm.pdb
	awk '$$1 == "ATOM" && $$9 > -15.5 && $$9 < 15.5 \
	    {a[$$6]++}END{for(pos in a) print pos}' $< > $@

# calculate SASA
ancillary/3sn6-sasa.tsv: ancillary/pdb/3sn6.pdb
	freesasa --format=seq $< | \
	    awk '$$2 == "R" && $$3 < 1000{print $$3,$$6}' OFS='\t' > $@
