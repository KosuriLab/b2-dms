"""
Simple wrapper to call consensus from Jensen-Shannon output
Nathan Lubock
"""

# Ensure Python 2/3 compatibility
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import argparse
import sys
import csv
from signal import signal, SIGPIPE, SIG_DFL
from collections import Counter

# catch broken pipe errors to allow ex) python pyParse.py foo bar | head
# see: https://stackoverflow.com/a/30091579
signal(SIGPIPE, SIG_DFL)

#===============================================================================

def consensus(string):
    """
    Generate a consensus sequence from an input string. If there are mo

    Parameters:
    -----------
    string :: str
        String representing all of the basecalls at that position

    Returns:
    --------
    'X' - tie at this position
    Otherwise the consensus sequence

    Requires:
    ---------
    Counter
    """
    uniqd = Counter(string).most_common()
    if len(uniqd) > 1 and uniqd[0][1] == uniqd[1][1]:
        return('X')
    else:
        return(uniqd[0][0])

def restricted_float(x):
    """
    Simple function to ensure argparse floats are within 0 - 1
    https://stackoverflow.com/a/12117065
    """
    x = float(x)
    if x < 0.0 or x > 1.0:
        raise argparse.ArgumentTypeError("{} not in range [0.0, 1.0]".format(x))
    return(x)

#===============================================================================

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Read Jensen-Shannon output and generate consensus sequence')
    parser.add_argument('infile',
                        type=argparse.FileType('r'),
                        default=sys.stdin,
                        nargs='?',
                        help='path to Jensen-Shannon file (or stdin if none)')
    parser.add_argument('-t',
                        '--gap-thresh',
                        dest='gap_thresh',
                        type=restricted_float,
                        default=0.5,
                        metavar='N',
                        help='ignore gaps if they are < this percent of sequence')
    args = parser.parse_args()

    #---------------------------------------------------------------------------

    tsv = csv.reader(args.infile, delimiter = str(u'\t'))
    for row_num, score, align in tsv:
        if float(align.count('-')) / len(align) < args.gap_thresh:
            cons = consensus(''.join(x for x in align if x != '-'))
        else:
            cons = consensus(align)
        print('{}\t{}'.format(cons, score), file=sys.stdout)


