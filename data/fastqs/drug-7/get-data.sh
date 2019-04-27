#!/bin/bash
for f in /data/projects/GPCR/raw-data/02162018_Drug_More-RNA_NextSeq/*.fastq.gz; do ln -s $f $(basename $f); done
mv 1_S1_R1_001.fastq.gz 0_1.fastq.gz
mv 2_S2_R1_001.fastq.gz 0_2.fastq.gz
mv 3_S3_R1_001.fastq.gz 150_1.fastq.gz
mv 4_S4_R1_001.fastq.gz 150_2.fastq.gz
mv 5_S5_R1_001.fastq.gz 625_1.fastq.gz
mv 6_S6_R1_001.fastq.gz 625_2.fastq.gz
mv 7_S7_R1_001.fastq.gz F_1.fastq.gz
mv 8_S8_R1_001.fastq.gz F_2.fastq.gz
touch rename.err
