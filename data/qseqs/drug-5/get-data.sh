#!/bin/bash
for f in ../../../../raw-data/11012017_Drug_ReReun/*_qseq.txt.gz; do ln -s $f $(basename $f); done
rename -v 's/s_1_([1-2])_([0-9]*).*/lane-1_$2_$1.qseq.gz/' *.gz > rename.err
