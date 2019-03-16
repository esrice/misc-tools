#!/usr/bin/env python3

import sys
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='Use the output of samtools '
            'mpileup to calculate pi and other statistics in sliding windows.')
    parser.add_argument('-c', '--min-coverage', type=int, default=10,
            help='minimum read coverage required to consider a position [10]')
    parser.add_argument('-C', '--max-coverage', type=int, default=100,
            help='maximum read coverage required to consider a position [100]')
    parser.add_argument('-w', '--window-size', type=int, default=20000,
            help='window size, in base pairs [20000]')
    # TODO empirically choose good defaults for these cutoffs
    parser.add_argument('-m', '--min-het-pct', type=float, default=0.25,
            help='minimum portion of reads supporting alternate allele for '
            'a position to be considered heterozygous [0.25]')
    parser.add_argument('-M', '--max-het-pct', type=float, default=0.75,
            help='maximum portion of reads supporting alternate allele for '
            'a position to be considered heterozygous [0.75]')
    return parser.parse_args()

# TODO make sure this function actually does what I want
# it to, because I copypastad it from a script used for a
# different purpose
def count_refs(pileup_base_string):
    ref_count = 0
    i = 0
    while i < len(pileup_base_string):
        if pileup_base_string[i] in ['.', ',']:
            ref_count += 1
        # if the current character is '^', then the next
        # character will be a PHRED quality, so skip it.
        elif pileup_base_string[i] in ['^']:
            i += 1
        i += 1

        # for other things, like insertions, deletions, and
        # ends of reads, there will not be any '.' or ','
        # characters to confuse us, so we don't need to
        # actually parse these out since all we care about
        # is the number of ref alleles

    return ref_count

def parse_pileup(instream, min_coverage, max_coverage):
    """
    Given a stream of the output of samtools mpileup, count
    the number of reads at each position supporting the
    reference allele and the number supporting a
    non-reference allele.

    Arguments:
    - instream: input file stream containing samtools
      mpileup output
    - min_coverage, max_coverage: the minimum and maximum
      coverage (inclusive) at a position to consider it

    Yields: tuples containing:
    - chromosome/sequence name
    - position in sequence
    - number of reads supporting reference
    - number of reads not supporting reference
    """
    for line in instream:
        splits = line.strip().split('\t')
        coverage = int(splits[3])

        if coverage < min_coverage or coverage > max_coverage:
            continue

        chrom = splits[0]
        position = int(splits[1])

        number_ref_reads = count_refs(splits[4])
        number_alt_reads = coverage - number_ref_reads

        yield (chrom, position, number_ref_reads, number_alt_reads)

def main():
    args = parse_args()

    for pile in parse_pileup(sys.stdin, args.min_coverage, args.max_coverage):
        chrom, position, number_ref_reads, number_alt_reads = pile
        print('\t'.join(map(str, [chrom, position, number_ref_reads,
            number_alt_reads])))
    # TODO calculate stats for windows

if __name__ == "__main__":
    main()
