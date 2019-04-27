"""
Simple barcode mapping utility.
Nathan Lubock
"""

# Ensure Python 2/3 compatibility
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import itertools
import argparse
import sys
from signal import signal, SIGPIPE, SIG_DFL

# catch broken pipe errors to allow ex) python pyParse.py foo bar | head
# see: https://stackoverflow.com/a/30091579
signal(SIGPIPE, SIG_DFL)

#===============================================================================

def fasta_reader(fasta):
    """
    Read in a fasta file lazily and return a generator of the name and sequence

    Parameters:
    -----------
    fasta :: FileType
        opened file

    Yields:
    -------
    generator :: (name, seq)
        name :: str
            Name of the read taken from the fasta file
        read :: str
            Sequence taken from the fasta file

    Requires:
    ---------
    itertools

    Example:
    --------
    itertools.groupby takes a key function and groups all items into a list
    until that key changes. We can key on lines beginning with >, then grab
    every line until the next record in the fasta. This makes our method robust
    to some fasta formats that have forced line breaks at given characters.

    foo = '>ABC>DEF>GHI'
    [(k, list(g)) for k,g in itertools.groupby(foo, lambda x: x == '>')]
    --> [(True, ['>']), (False, ['A', 'B', 'C']), (True, ['>']), ... ]

    Note:
    -----
    Adapted from: https://www.biostars.org/p/710/#1412
    """
    # ditch the boolean (x[0]) and just keep the header/seq grouping
    fa_iter = (x[1] for x in itertools.groupby(fasta, lambda line: line[0] == ">"))
    for header in fa_iter:
        # drop the ">"
        name = next(header)[1:].strip()
        # join all sequence lines to one by iterating until the next group.
        read = "".join(s.strip() for s in next(fa_iter))
        yield name, read

#-------------------------------------------------------------------------------

def revcomp(seq):
    """
    Reverse Complement a string

    Parameters:
    -----------
    seq :: str

    Returns:
    --------
    str
    """
    comp = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N'}
    return ''.join(comp[nuc] for nuc in seq[::-1])

#===============================================================================

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Generate all frameshifts for a DNA fasta file')
    parser.add_argument('infile',
                        type=argparse.FileType('r'),
                        default=sys.stdin,
                        nargs='?',
                        help='path to a *.fasta file of the reads (or stdin if none)')
    parser.add_argument('-r',
                        '--rev-comp',
                        dest='rc',
                        action='store_true',
                        help='reverse complement the sequence?')
    args = parser.parse_args()

    for head, seq in fasta_reader(args.infile):
        if args.rc:
            seq = revcomp(seq)
        out = ((head + '_' + str(x + 1), seq[:x] + seq[x + 1:])
                for x in xrange(0, len(seq)))
        for header, dels in out:
            print('>{}\n{}'.format(header, dels), file=sys.stdout)




