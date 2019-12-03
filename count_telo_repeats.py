#!/usr/bin/env python3
"""
Count the number of telomeric repeats in the final bases of each
sequence in a fasta file, and output a table containing sequence name,
telomeric repeat count, and p-value.
"""

import argparse
import re

from scipy.stats import binom
import screed


telomeric_repeat = re.compile('TTAGGG', re.IGNORECASE)

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-b', '--bases', type=int, default='1000',
                        help='number of bases at end of sequence to search')
    parser.add_argument('fasta', type=screed.open, help='fasta file '
                        'containing sequences to search for telomeres')
    return parser.parse_args()


def count_telo_repeats(seq_to_count):
    """
    Count number of occurrences of 'TTAGGG' in a sequence.
    """
    return len(telomeric_repeat.findall(seq_to_count))


def main():
    args = parse_args()
    p_TTAGGG = 0.25**6

    print('\t'.join(['seq_name', 'telo_repeat_count', 'p_val']))
    for seq in args.fasta:
        num_repeats = count_telo_repeats(seq.sequence[-1*args.bases:])
        p = 1.0 - binom.cdf(num_repeats, args.bases, p_TTAGGG)
        print('\t'.join(map(str, [seq.name, num_repeats, p])))


if __name__ == '__main__':
    main()
