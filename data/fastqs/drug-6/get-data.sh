#!/bin/bash
for f in /data/projects/GPCR/raw-data/11152017_Drug_NextSeq/*.fastq.gz; do ln -s $f $(basename $f); done
rename 's/_S._L00(.)_R1_001/_$1/' *.gz > rename.err
