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
import multiprocessing
import sys
import csv
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

#-------------------------------------------------------------------------------

def translate(pack):
    """
    Translate a DNA sequence into it's protein equivalent

    Parameters:
    -----------
    pack :: (str, str, bool)
        var - variant header
        seq - DNA sequnce
        rc - reverse complement the sequence before translating?

    Returns:
    --------
    (var, aminos, flag) :: (str, str, bool)
        var - variant header
        aminos - Protein sequence

    Requires:
    ---------
    itertools

    Note:
    -----
    Assumes sequences are already in-frame (e.g. does not find start codons).
    Will convert any N's to a dummy codon.
    Adapted from: http://stackoverflow.com/a/19522901
    """
    codon_table = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W'}
    var, seq, rc = pack
    if rc:
        dna = revcomp(seq)
    else:
        dna = seq

    codons = (dna[x:x+3] for x in range(0, len(dna), 3))
    coding_seq = itertools.takewhile(lambda x: len(x) == 3, codons)
    clean_codons = ('X' if 'N' in x else codon_table[x] for x in coding_seq)
    aminos = ''.join(clean_codons)

    # set a flag for downstream filtering
    if 'X' in aminos:
        return (var, aminos, True)
    else:
        return (var, aminos, False)

#===============================================================================

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Translate DNA fasta to Protein. Converts N to X and STOPs to *')
    parser.add_argument('infile',
                        type=argparse.FileType('r'),
                        default=sys.stdin,
                        nargs='?',
                        help='path to a *.fasta file of the reads (or stdin if none)')
    parser.add_argument('-v',
                        '--verbose',
                        dest='verbose',
                        action='store_true',
                        help='Verbose logging')
    parser.add_argument('-j',
                        '--proc',
                        dest='proc',
                        type=int,
                        default=1,
                        metavar='N',
                        choices=range(1, multiprocessing.cpu_count()+1),
                        help='number of processors (default=1, max={})'.format(
                            multiprocessing.cpu_count()))
    parser.add_argument('-b',
                        '--bad-bcs',
                        dest='bad_bcs',
                        type=argparse.FileType('w'),
                        metavar='foo',
                        help='output bad barcodes to this file')
    parser.add_argument('-r',
                        '--rev-comp',
                        dest='rc',
                        action='store_true',
                        help='reverse complement the sequence before translating?')
    args = parser.parse_args()

    #---------------------------------------------------------------------------

    if args.verbose:
        print('Translating...', file=sys.stderr)

    pool = multiprocessing.Pool(args.proc)
    pack = ((head, seq, args.rc) for head, seq in fasta_reader(args.infile))
    raw_prots = pool.imap_unordered(translate, pack, chunksize=10000)
    # set up writers
    bad_flag = False
    if args.bad_bcs is not None:
        bad_out = csv.writer(args.bad_bcs, lineterminator='\n')
        bad_flag = True

    bad_count = 0
    good_count = 0
    for header, prot, bad_prot in raw_prots:
        if bad_prot:
            bad_count += 1
            if bad_flag:
                bad_out.writerow((header, prot))
        else:
            print(">{}\n{}".format(header, prot), file=sys.stdout)
            good_count += 1

    if args.verbose:
        print('Removed {} BCs with a N base call in variant'.format(bad_count),
              file=sys.stderr)
        print('Kept {} BCs'.format(good_count), file=sys.stderr)

    pool.close()
    pool.join()
