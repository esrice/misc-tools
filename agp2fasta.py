#!/usr/bin/env python3
"""
Given contigs in fasta format and their order and orientation into
scaffolds in AGP format, outputs the assembled scaffolds in fasta
format.
"""

import argparse
import sys

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna

import agp

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-o', '--outfile', type=argparse.FileType('w'),
                        help='where to write fasta of scaffolds [stdout]',
                        default=sys.stdout)
    parser.add_argument('contigs_fasta',
                        type=lambda f: SeqIO.parse(f, 'fasta'),
                        help='Assembly to flip scaffolds in')
    parser.add_argument('agp', type=lambda s: agp.read_agp(open(s)),
                        help='AGP file assembling contigs into scaffolds')
    return parser.parse_args()


def main():
    args = parse_args()

    # unfortunately, the contigs fasta file I'm writing this for has
    # variable line-length and is thus not faidx-able, so we have to
    # load it into memory :(
    contigs = {record.id: record.seq for record in args.contigs_fasta}

    current_sequence = None
    current_chrom = None
    # loop through AGP, skipping comment lines, which my agp library
    # yields as strings
    for row in filter(lambda r: not isinstance(r, str), args.agp):
        # check if starting a new chromosome
        if row.object != current_chrom:
            # if this is not the first chromosome, output the previous
            # chromosome
            if current_chrom is not None:
                record = SeqRecord(current_sequence, id=current_chrom,
                                   description='')
                print(record.format('fasta'), end='')
            # start the new chromosome as an empty sequence
            current_chrom = row.object
            current_sequence = Seq('', generic_dna)

        if row.is_gap:
            current_sequence += 'N' * row.gap_length
        else:
            start, end = row.component_beg - 1, row.component_end
            component = contigs[row.component_id][start:end]
            if row.orientation == '-':
                component = component.reverse_complement()
            current_sequence += component

    record = SeqRecord(current_sequence, id=current_chrom, description='')
    print(record.format('fasta'), end='')


if __name__ == '__main__':
    main()
