#!/usr/bin/env python3
"""
Given an AGP file describing an assembly and a list of chromosomes to
flip, output a new AGP file with those chromosomes reoriented in
reverse-complement.
"""

import argparse
import sys

# names of fields in an AGP file, in order, for later use as dictionary
# keys
AGP_FIELD_NAMES = [
    'object',
    'object_beg',
    'object_end',
    'part_number',
    'component_type',
    'component_id',
    'component_beg',
    'component_end',
    'orientation',
]

class AgpRow:
    """
    A single non-comment row of an AGP file. Because AGP is a weird
    kind of TSV where fields can have different meanings depending on
    the value of the component_type field, I'm representing this as a
    class with instance variables rather than an OrderedDict or
    something along those lines. See the NCBI documentation for more
    information:

    https://www.ncbi.nlm.nih.gov/assembly/agp/AGP_Specification/
    """
    def __init__(self, line):
        """
        Creates a new instance of AgpRow by parsing the given line of
        text.
        """
        splits = line.strip().split('\t')
        self.object = splits[0]
        self.object_beg = int(splits[1])
        self.object_end = int(splits[2])
        self.part_number = int(splits[3])
        self.component_type = splits[4]

        if self.component_type in ['N', 'U']:
            self.is_gap = True
            self.gap_length = int(splits[5])
            self.gap_type = splits[6]
            self.linkage = splits[7]
            self.linkage_evidence = splits[8]
        else:
            self.is_gap = False
            self.component_id = splits[5]
            self.component_beg = int(splits[6])
            self.component_end = int(splits[7])
            self.orientation = splits[8]

    def __str__(self):
        """
        Returns the string representation of the AGP row as a line of
        text containing all the fields separated by tabs.
        """
        if self.is_gap:
            return '\t'.join(map(str, [self.object, self.object_beg,
                                       self.object_end, self.part_number,
                                       self.component_type, self.gap_length,
                                       self.gap_type, self.linkage,
                                       self.linkage_evidence]))
        else:
            return '\t'.join(map(str, [self.object, self.object_beg,
                                       self.object_end, self.part_number,
                                       self.component_type, self.component_id,
                                       self.component_beg, self.component_end,
                                       self.orientation]))


def read_agp(infile):
    """
    Reads an AGP file, yielding rows as OrderedDict objects with field
    key names specified by AGP_FIELD_NAMES. Yields comment lines as
    strings.
    """
    for line in infile:
        if line.startswith('#'):
            yield line.strip()
        else:
            yield AgpRow(line)


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
    parser.add_argument('agp', nargs='?', type=lambda s: read_agp(open(s)),
                        help='AGP file to modify [STDIN]',
                        default=read_agp(sys.stdin))
    return parser.parse_args()


def reverse_agp(chromosome_rows):
    """
    Given a list of AGP rows assembling a complete chromosome, return a
    new set of rows assembling the reverse complement of the input rows.

    Args:
        chromosome_rows (list): a list of OrderedDict objects, each
            containing one AGP entry with keys specified by
            AGP_FIELD_NAMES

    Returns:
        reverse_chromosome_rows (list): the input, except in reverse
            order, with the rows modified to be reverse-complements
    """
    # first, get the total number of elements in the chromosome and the
    # total length of the chromosome so we can count backwards from
    # there.
    number_parts = chromosome_rows[-1].part_number
    chrom_length = chromosome_rows[-1].object_end

    reversed_rows = []
    for row in reversed(chromosome_rows):
        # TODO debug this carefully to look for OBOEs
        new_start = chrom_length - row.object_end + 1
        new_end = chrom_length - row.object_beg + 1
        row.object_beg, row.object_end = new_start, new_end

        row.part_number = number_parts - row.part_number + 1

        if not row.is_gap:
            if row.orientation == '+':
                row.orientation = '-'
            elif row.orientation == '-':
                row.orientation = '+'
            # if the orientation is something else, leave it be

        reversed_rows.append(row)

    return reversed_rows

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
                    rows_to_print = reverse_agp(rows_this_chromosome)
                else:
                    rows_to_print = rows_this_chromosome

                for row_to_print in rows_to_print:
                    print(row_to_print, file=args.outfile)

                rows_this_chromosome = []

            rows_this_chromosome.append(row)

    # finish up the last chromosome
    if rows_this_chromosome[-1].object in args.chroms_to_flip:
        rows_to_print = reverse_agp(rows_this_chromosome)
    else:
        rows_to_print = rows_this_chromosome

    for row_to_print in rows_to_print:
        print(row_to_print, file=args.outfile)

if __name__ == '__main__':
    main()
