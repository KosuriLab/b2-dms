for f in ../../../../raw-data/070717_InsLib_ForskSanityCheck/*.txt.gz; do ln -s $f $(basename $f); done
rename -v 's/s_[0-9]_([1-2])_([0-9]*).*/rep-1_$2_$1.qseq.gz/' *.gz > rename.err

