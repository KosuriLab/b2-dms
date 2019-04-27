"""
dssp.py - quick wrapper to process dssp files
Nathan Lubock
"""

import itertools
import argparse
import multiprocessing
import sys
import csv
from signal import signal, SIGPIPE, SIG_DFL

# EXTERNAL DEPS
# Note dms_tools requires at least py3.4
# If you have disulfide bonds in your protein, you must use commit:
# 34228b9556150ed8a9d10c57ddcbee115022b6c2 or later
from dms_tools2 import dssp

# catch broken pipe errors
# see: https://stackoverflow.com/a/30091579
signal(SIGPIPE, SIG_DFL)

#===============================================================================

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Parse DSSP files for relevant information')
    parser.add_argument('infile',
                        type=str,
                        help='path to a *.dssp file')
    parser.add_argument('-c',
                        '--chain',
                        dest='chain',
                        nargs='?',
                        help='PDB chain to process (default all)')
    args = parser.parse_args()

    df = dssp.processDSSP(args.infile, args.chain)
    df.sort_values(by='site').to_csv(sys.stdout, sep='\t', index=False)


