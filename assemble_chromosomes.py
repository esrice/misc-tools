#!/usr/bin/env python3

import argparse
import sys

def parse_args():
    parser = argparse.ArgumentParser(description='Given a fasta containing '
            'scaffolds and a file laying out scaffolds into chromosomes, '
            'output a fasta of chromosomes.')
    parser.add_argument('scaffolds_fasta', help='a fasta file containing '
            'scaffolds to be assembled into chromosomes')
    parser.add_argument('layout_file', type=argparse.FileType('r'),
            help='a file with two columns: first is the name of a chromosome, '
            'second is a comma-separated list of scaffolds in that chromosome '
            'in the format "[name]:[orientation]", e.g. "scaf5:+,scaf2:-" '
            'indicates that this chromosome is made up of scaf5 first in the '
            'forward orientation, followed by scaf2 in the reverse complement '
            'orientation.')
    parser.add_argument('-p', '--path-to-samtools', default='samtools',
            help='the path to the samtools executable, if it\'s not in $PATH')
    return parser.parse_args()

"""
Named tuple representing a scaffold and its orientation.
Fields:
    - name: ID of scaffold as it appears in fasta header
    - orientation: true if (+), false if (-)
"""
Scaffold = namedtuple('Scaffold', ['name', 'orientation'])

def parse_scaffold_list(scaffold_list_string):
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
    for scaffold in scaffold_strings:
        name, orientation = scaffolds.split(':')
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
    used_scaffolds = set()

    for line_number, line in enumerate(args.layout_file):
        try:
            chrom, scaffolds = line.strip().split()
            scaffolds = parse_scaffolds_list(scaffolds)
        except ValueError:
            print('FATAL: Misformatted layout file on line {}'
                    .format(line_number + 1), file=sys.stderr)
            sys.exit(1)

        chromosome_sequence = ''
        for scaffold_number, scaffold in enumerate(scaffolds):
            # TODO start here
            # 1. Append some 'N's if scaffold_number != 0
            # 2. Get sequence using samtools faidx
            # 3. Append sequence to chromosome
        # 4. format and print fasta entry

if __name__ == '__main__':
    main()
