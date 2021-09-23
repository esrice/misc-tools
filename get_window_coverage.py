#!/usr/bin/env python3

"""
Use the output of samtools mpileup to calculate the average coverage
per window.

Example full usage:
samtools mpileup -f ref.fa -q30 aligned.bam | get_window_coverage.py > windows.tsv
"""

import argparse
import sys


def parse_args():
    """Parse command-line arguments"""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-w",
        "--window-size",
        type=int,
        default=50000,
        help="window size, in base pairs [50000]",
    )
    return parser.parse_args()


def parse_pileup(instream):
    """
    Given a stream of the output of samtools mpileup, count
    the number of reads at each position supporting the
    reference allele and the number supporting a
    non-reference allele.

    Arguments:
    - instream: input file stream containing samtools
      mpileup output

    Yields: tuples containing:
    - chromosome/sequence name
    - position in sequence
    - coverage at this position
    """
    for line in instream:
        splits = line.strip().split("\t")
        coverage = int(splits[3])
        chrom = splits[0]
        position = int(splits[1])

        yield (chrom, position, coverage)


def print_window(chromosome, window_start, window_end, sum_coverage):
    """
    Calculate the mean coverage for this window and print a bed-ish
    line containing the window boundaries and this mean coverage.

    Args:
        chromosome: name of chromosome (str) in which window resides
        window_start: beginning position of window (int)
        window_end: end position of window (int)
        sum_coverage: total number of bases mapping to this window

    Returns: nothing, but prints a tab-separated line to stdout with columns:
        * chromosome name
        * window start
        * window end position
        * mean coverage in window
    """
    print(
        "\t".join(
            map(
                str,
                [
                    chromosome,
                    window_start,
                    window_end,
                    sum_coverage / (window_end - window_start),
                ],
            )
        )
    )


def main():
    """main method of script"""
    args = parse_args()

    window_start = 1  # coordinate for the start of the current window
    sum_coverage = 0  # running sum of coverage for each base in this window
    last_chrom = ""  # name of last chromosome, to see if we're on a new one
    last_position = 0  # last position encountered, for use when switching chromosomes
    for chrom, position, coverage in parse_pileup(sys.stdin):
        # we are on a new chromosome, so report stats for the last window on
        # the previous chromosome
        if chrom != last_chrom:
            if last_chrom != "":
                print_window(last_chrom, window_start, last_position, sum_coverage)
            last_chrom = chrom
            window_start = position
            sum_coverage = 0

        # we're outside of the range of the last window, so report stats for
        # previous window and start a new one
        if position > window_start + args.window_size:
            print_window(chrom, window_start, position - 1, sum_coverage)
            window_start = position
            sum_coverage = 0

        sum_coverage += coverage
        last_position = position

    # don't forget the last window!
    print_window(last_chrom, window_start, last_position, sum_coverage)


if __name__ == "__main__":
    main()
