"""
synon-filter.py - simple script to check and trim reads
Nathan Lubock

Only outupt  barcodes that aligns exactly to our chunk. For a few chunks, trim
to the first whole codon if the nucleotides are correct.
"""

# Ensure Python 2/3 compatibility
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import argparse
import sys
from signal import signal, SIGPIPE, SIG_DFL

# catch broken pipe errors to allow ex) python pyParse.py foo bar | head
# see: https://stackoverflow.com/a/30091579
signal(SIGPIPE, SIG_DFL)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Filter .sam file before translating')
    parser.add_argument('sam',
                        type=argparse.FileType('r'),
                        default=sys.stdin,
                        nargs='?',
                        help='path to sam file (or stdin if none)')
    args = parser.parse_args()

    good_loc = set((37, 190, 354, 519, 678, 817, 982, 1132))
    trim = {354:'C', 519:'G', 678:'C'}

    for line in args.sam:
        split = line.split('\t')
        bc = split[0]
        pos = int(split[3])
        seq = split[9]
        if pos in good_loc:
            if pos in trim:
                if seq[0] == trim[pos]:
                    print('>{}\n{}'.format(bc, seq[1:]), file=sys.stdout)
            else:
                print('>{}\n{}'.format(bc, seq), file=sys.stdout)





