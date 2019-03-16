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
    parser.add_argument('-w', '--window-size', type=int, default=50000,
            help='window size, in base pairs [50000]')
    parser.add_argument('-p', '--min-fraction-window', type=float, default=0.2,
            help='minimum fraction of a window with sufficient coverage '
            'required in order to report window [0.2]')
    parser.add_argument('-m', '--min-het-fct', type=float, default=0.2,
            help='minimum fraction of reads supporting alternate allele for '
            'a position to be considered heterozygous [0.2]')
    parser.add_argument('-M', '--max-het-fct', type=float, default=0.8,
            help='maximum fraction of reads supporting alternate allele for '
            'a position to be considered heterozygous [0.8]')
    return parser.parse_args()

def count_refs(pileup_base_string):
    """
    Given a pileup string from samtools mpileup output,
    count the number of reads that support the reference
    allele in that position.
    """
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

def print_window_stats(hom_ref_count, het_count, hom_alt_count, chrom,
        window_start, window_end, min_fraction_window):
    """
    Calculate rates of heterozygous and homozygous
    non-reference positions in a window and print these
    to STDOUT in bed-ish format.

    Arguments:
    - hom_ref_count, het_count, hom_alt_count: number of
      positions in this window that are homozygous for
      the reference allele, heterozygous, or homozygous for
      an alternate allele
    - window_start, window_end: start/end position of
      window, in base pairs
    - args: args from parse_args function for getting
      parameters like window size

    Prints: tab-separated line with fields: chromosome
        name, window start position, window end position,
        rate of heterozygosity in window, rate of
        homozygous non-reference positions in window.
    """
    num_called_positions = hom_ref_count + het_count + hom_alt_count
    window_size = window_end - window_start
    if num_called_positions/window_size > min_fraction_window:
        het_rate = het_count/num_called_positions
        hom_alt_rate = hom_alt_count/num_called_positions
        print('\t'.join(map(str, [chrom, window_start, window_end,
            het_rate, hom_alt_rate])))

def main():
    args = parse_args()

    window_start = 1 # coordinate for the start of the current window
    # running counts of number of times we observe positions matching these
    # descriptions:
    hom_ref_count, het_count, hom_alt_count = 0, 0, 0
    for pile in parse_pileup(sys.stdin, args.min_coverage, args.max_coverage):
        chrom, position, number_ref_reads, number_alt_reads = pile

        # we're outside of the range of the last window, so report stats for
        # previous window and start a new one
        if position > window_start + args.window_size:
            print_window_stats(hom_ref_count, het_count, hom_alt_count, chrom,
                    window_start, position-1, args.min_fraction_window)
            window_start = position
            hom_ref_count, het_count, hom_alt_count = 0, 0, 0

        # calculate the fraction of reads supporting the alternate allele and
        # count this position accordingly
        alt_fraction = number_alt_reads/(number_ref_reads + number_alt_reads)
        if alt_fraction < args.min_het_fct:
            hom_ref_count += 1
        elif alt_fraction < args.max_het_fct:
            het_count += 1
        else:
            hom_alt_count += 1

    # don't forget the last window!
    print_window_stats(hom_ref_count, het_count, hom_alt_count, chrom,
            window_start, position, args.min_fraction_window)

if __name__ == "__main__":
    main()
