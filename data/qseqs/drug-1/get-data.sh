#!/bin/bash
for f in ../../../../raw-data/061117_HiSeq/*_qseq.txt.gz; do ln -s $f $(basename $f); done
rename -v 's/s_1_([1-2])_([0-9]*).*/rep-1_$2_$1.qseq.gz/' *.gz > rename.err
for f in ../../../../raw-data/061217_InsLibMM1_pilot/*_qseq.txt.gz; do ln -s $f $(basename $f); done
rename -v 's/s_1_([1-2])_([0-9]*).*/rep-2_$2_$1.qseq.gz/' *.gz >> rename.err
