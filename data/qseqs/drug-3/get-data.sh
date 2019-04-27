#!/bin/bash
for f in /data/projects/GPCR/raw-data/092017_Forsk-Iso_Lane-1/*_qseq.txt.gz; do ln -s $f $(basename $f); done
rename -v 's/s_1_([1-2])_([0-9]*).*/lane-1_$2_$1.qseq.gz/' *.gz > rename.err
for f in /data/projects/GPCR/raw-data/092017_Forsk-Iso_Lane-2/*_qseq.txt.gz; do ln -s $f $(basename $f); done
rename -v 's/s_2_([1-2])_([0-9]*).*/lane-2_$2_$1.qseq.gz/' *.gz >> rename.err
