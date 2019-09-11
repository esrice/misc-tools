#!/usr/bin/env python3

import argparse
from collections import defaultdict
import textwrap

import screed


def breaks_type(filename):
    with open(filename, 'r') as breaks_file:
        contigs_to_break = defaultdict(list)
        for line in (l for l in breaks_file if not l.startswith('#')):
            splits = line.strip().split()
            contigs_to_break[splits[0]].append((splits[1],
                                                int(splits[2]),
                                                int(splits[3])))
        return contigs_to_break


def format_fasta(name, seq, wrap=60):
    """
    Formats a sequence as a string of a fasta entry.

    >>> print(format_fasta('seq_name', 'ACTGTCATCATGTGC', wrap=5))
    >seq_name
    ACTGT
    CATCA
    TGTGC
    """
    return ">{}\n{}".format(name, textwrap.fill(seq, width=wrap))


def break_seq(seq, breaks):
    """
    Breaks the sequence into pieces with names and boundaries specified
    by breaks.

    Args:
        seq (str): String to break into pieces
        breaks ((name (str), start (int), end (int))): name, start
            position, and end position for each chunk

    Returns:
        Generator over (break_name (str), break_seq (str)) pairs

    >>> seq = 'ACTAGCTAGTCGACAATTATCGCGACT'
    >>> breaks = [('A', 0, 5), ('B', 5, 17), ('C', 17, 27)]
    >>> list(break_seq(seq, breaks))
    [('A', 'ACTAG'), ('B', 'CTAGTCGACAAT'), ('C', 'TATCGCGACT')]
    """
    for break_name, break_start, break_end in breaks:
        yield (break_name, seq[break_start:break_end])


def parse_args():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('breaks', type=breaks_type, help='tsv containing '
                        'the following fields: contig name, subcontig name, '
                        'subcontig start, subcontig end')
    parser.add_argument('fasta', type=screed.open, help='fasta file '
                        'containing contigs to break')
    return parser.parse_args()


def main():
    args = parse_args()

    for seq in args.fasta:
        if seq.name in args.breaks:
            new_seqs = break_seq(seq.sequence, args.breaks[seq.name])
            for new_seq_name, new_seq_sequence in new_seqs:
                print(format_fasta(new_seq_name, new_seq_sequence))
        else:
            print(format_fasta(seq.name, seq.sequence))


if __name__ == '__main__':
    main()
