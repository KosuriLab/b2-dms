"""
Generate all synonymous mutants for a input protein
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

#-------------------------------------------------------------------------------

def to_prot(dna):
    """
    Translate a DNA sequence into it's protein equivalent

    Parameters:
    -----------
    dna - DNA sequnce

    Returns:
    --------
    (aminos, flag) :: (str, bool)
        aminos - Protein sequence
        flag - Was there an X?

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
    codons = (dna[x:x+3] for x in range(0, len(dna), 3))
    coding_seq = itertools.takewhile(lambda x: len(x) == 3, codons)
    clean_codons = ('X' if 'N' in x else codon_table[x] for x in coding_seq)
    aminos = ''.join(clean_codons)

    # set a flag for downstream filtering
    if 'X' in aminos:
        return (aminos, True)
    else:
        return (aminos, False)

#-------------------------------------------------------------------------------

def to_dna(dna):
    rev_codon_table = {
    'A':['GCT','GCC','GCA','GCG'],
    'R':['CGT','CGC','CGA','CGG','AGA','AGG'],
    'N':['AAT','AAC'],
    'D':['GAT','GAC'],
    'C':['TGT','TGC'],
    'Q':['CAA','CAG'],
    'E':['GAA','GAG'],
    'G':['GGT','GGC','GGA','GGG'],
    'H':['CAT','CAC'],
    'I':['ATT','ATC','ATA'],
    'L':['TTA','TTG','CTT','CTC','CTA','CTG'],
    'K':['AAA','AAG'],
    'M':['ATG'],
    'F':['TTT','TTC'],
    'P':['CCT','CCC','CCA','CCG'],
    'S':['TCT','TCC','TCA','TCG','AGT','AGC'],
    'T':['ACT','ACC','ACA','ACG'],
    'W':['TGG'],
    'Y':['TAT','TAC'],
    'V':['GTT','GTC','GTA','GTG'],
    '*':['TAA','TGA','TAG']}
    # get codons
    old_codons = (dna[x:x+3] for x in range(0, len(dna), 3))
    coding_seq = list(itertools.takewhile(lambda x: len(x) == 3, old_codons))
    prot, flag = to_prot(dna)

    if flag:
        raise ValueError('N found in input sequence')

    raw_seqs = list()
    names = list()
    for pos, aa in enumerate(prot):
        new_codons = rev_codon_table[aa]
        for codon in new_codons:
            names.append('{}_{}_{}'.format(pos + 1, aa, codon))
            raw_seqs.append(''.join(
                itertools.chain.from_iterable([coding_seq[:pos], codon, coding_seq[pos+1:]])))

    return([x for x in zip(names, raw_seqs) if x[1] != dna])



#===============================================================================

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Generate all single synonymous mutations for an input sequence')
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

    # drop the fasta header since we don't need it
    for _, seq in fasta_reader(args.infile):
        if args.rc:
            seq = revcomp(seq)
        out = to_dna(seq)
        for header, synon in out:
            print('>{}\n{}'.format(header, synon), file=sys.stdout)


