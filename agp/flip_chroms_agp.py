#!/usr/bin/env python3
"""
Given an AGP file describing an assembly and a list of chromosomes to
flip, output a new AGP file with those chromosomes reoriented in
reverse-complement.
"""

import argparse
import sys

import agp

def chroms_list_type(arg):
    """
    argparse type function that opens the file at the specified location
    and then returns a set of the stripped lines
    """
    return set(map(lambda l: l.strip(), open(arg, 'r')))

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-o', '--outfile', type=argparse.FileType('w'),
                        nargs='?', help='where to write output AGP [STDOUT]',
                        default=sys.stdout)
    parser.add_argument('chroms_to_flip', type=chroms_list_type,
                        help='list of chromosomes to flip, one per line')
    parser.add_argument('agp', nargs='?', type=lambda s: agp.read_agp(open(s)),
                        help='AGP file to modify [STDIN]',
                        default=agp.read_agp(sys.stdin))
    return parser.parse_args()


def main():
    args = parse_args()

    rows_this_chromosome = []
    for row in args.agp:
        # comment lines are yielded by read_agp() as strings, so we
        # just print the line as is and move on
        if isinstance(row, str):
            print(row, file=args.outfile)
        else:
            if (rows_this_chromosome
                    and row.object != rows_this_chromosome[-1].object):
                # we are now on a new chromosome, so let's output all
                # rows for the previous one, reversing if necessary.
                if rows_this_chromosome[-1].object in args.chroms_to_flip:
                    rows_to_print = agp.reverse_agp(rows_this_chromosome)
                else:
                    rows_to_print = rows_this_chromosome

                for row_to_print in rows_to_print:
                    print(row_to_print, file=args.outfile)

                rows_this_chromosome = []

            rows_this_chromosome.append(row)

    # finish up the last chromosome
    if rows_this_chromosome[-1].object in args.chroms_to_flip:
        rows_to_print = agp.reverse_agp(rows_this_chromosome)
    else:
        rows_to_print = rows_this_chromosome

    for row_to_print in rows_to_print:
        print(row_to_print, file=args.outfile)


if __name__ == '__main__':
    main()
