#!/usr/bin/env python3

import argparse
import sys
import pyfaidx
from collections import namedtuple
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def parse_args():
    parser = argparse.ArgumentParser(description='Given a fasta containing '
            'scaffolds and a file laying out scaffolds into chromosomes, '
            'output a fasta of chromosomes.')
    parser.add_argument('scaffolds_fasta', type=pyfaidx.Fasta, help='a fasta '
            'file containing scaffolds to be assembled into chromosomes')
    parser.add_argument('layout_file', type=argparse.FileType('r'),
            help='a file with two columns: first is the name of a chromosome, '
            'second is a comma-separated list of scaffolds in that chromosome '
            'in the format "[name]:[orientation]", e.g. "scaf5:+,scaf2:-" '
            'indicates that this chromosome is made up of scaf5 first in the '
            'forward orientation, followed by scaf2 in the reverse complement '
            'orientation.')
    parser.add_argument('-g', '--gap-size', type=int, default=1000,
            help='number of N\'s to place between scaffolds')
    return parser.parse_args()

"""
Named tuple representing a scaffold and its orientation.
Fields:
    - name: ID of scaffold as it appears in fasta header
    - orientation: true if (+), false if (-)
"""
Scaffold = namedtuple('Scaffold', ['name', 'orientation'])

def parse_scaffolds_list(scaffold_list_string):
    """
    Given a comma-separated list of scaffolds and their
    orientations, turn it into a list of Scaffold objects

    scaffold_list_string: a string containing a comma-
        separated list of scaffolds and orientations in the
        format '[name]:[orientation]', e.g.,
        'scaf5:+,scaf2:-,scaf10:-'

    returns: a list of Scaffold objects
    """
    scaffold_strings = scaffold_list_string.split(',')
    scaffolds = []
    for scaffold_string in scaffold_strings:
        name, orientation = scaffold_string.split(':')
        if orientation == '+':
            orientation = True
        elif orientation == '-':
            orientation = False
        else:
            # the calling function will catch this and make
            # a nice error message with the line number
            raise ValueError

        scaffolds.append(Scaffold(name, orientation))

    return scaffolds

def main():
    args = parse_args()

    # set of string IDs of scaffolds that have been placed
    # onto chromosomes, so that all unplaced scaffolds can
    # be printed out at the end
    assigned_scaffolds = set()

    for line_number, line in enumerate(args.layout_file):
        try: # make sure there are the right number of columns
            chromosome_name, scaffolds = line.strip().split()
            scaffolds = parse_scaffolds_list(scaffolds)
        except ValueError:
            print('FATAL: Misformatted layout file on line {}'
                    .format(line_number + 1), file=sys.stderr)
            sys.exit(1)

        chromosome_sequence = ''
        for scaffold_number, scaffold in enumerate(scaffolds):
            # add a gap if this isn't the first scaffold
            if scaffold_number != 0:
                chromosome_sequence += 'N' * args.gap_size

            # get sequence using samtools faidx
            try:
                this_scaffold_sequence = args.scaffolds_fasta[scaffold.name]
            except KeyError:
                eprint('FATAL: no sequence with ID "{}" in input'.format(
                    scaffold.name))
                sys.exit(1)

            # add scaffold to chromosome
            if scaffold.orientation:
                chromosome_sequence += this_scaffold_sequence[:].seq
            else: # reverse complement sequence if necessary
                chromosome_sequence += (-this_scaffold_sequence[:]).seq

            # take a note that this scaffold has been
            # assigned to a chromosome
            assigned_scaffolds.add(scaffold.name)

        # format and print chromosome
        print(SeqRecord(Seq(chromosome_sequence), id=chromosome_name,
            description='').format('fasta'), end='')

    # loop through fasta file outputting all unplaced scaffolds
    for scaffold in args.scaffolds_fasta:
        if scaffold.name not in assigned_scaffolds:
            print(SeqRecord(Seq(scaffold[:].seq), id=scaffold.name,
                description='').format('fasta'), end='')

if __name__ == '__main__':
    main()
