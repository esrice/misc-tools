#!/usr/bin/env python3

"""
flip_chroms.py -- reverse-complement sequences in an
assembly.
"""

import argparse
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

def parse_args():
    parser = argparse.ArgumentParser(description='Reverse complement given '
            'sequences in an assembly.')
    parser.add_argument('assembly_fasta',
            type=lambda f: SeqIO.parse(f, 'fasta'),
            help='Assembly to flip scaffolds in')
    parser.add_argument('scaffolds_to_flip', type=argparse.FileType('r'),
            help='File containing seq IDs, one per line, to be revcomped')
    return parser.parse_args()

def main():
    args = parse_args()

    scaffolds_to_flip = set(map(lambda l: l.strip(), args.scaffolds_to_flip))

    for record in args.assembly_fasta:
        if record.id in scaffolds_to_flip:
            record.seq = record.seq.reverse_complement()

        print(record.format('fasta'), end='')

if __name__ == '__main__':
    main()
