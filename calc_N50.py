#!/usr/bin/env python3
"""
Given one or more fastq/fasta files, calculate the N50.
"""
import argparse
import sys

import screed

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('seqfiles', type=screed.open, nargs='+',
                        help='fasta/q[.gz] file(s) containing sequences')
    return parser.parse_args()

def main():
    args = parse_args()

    read_lengths = []
    for seq_file in args.seqfiles:
        for read in seq_file:
            read_lengths.append(len(read.sequence))

    total_length = sum(read_lengths)
    read_lengths.sort(reverse=True)

    cumulative_sum = 0
    for l in read_lengths:
        cumulative_sum += l
        if cumulative_sum >= total_length / 2:
            print(l)
            sys.exit()

if __name__ == '__main__':
    main()
