#!/usr/bin/env python3
"""
Given a sam/bam file, calculate the mean coverage for each sequence
therein. Print out a tsv where the three columns are reference
sequence name, mean alignment coverage over the length of that
sequence, and sequence length.
"""

import argparse

import pysam


def parse_args():
    """parse command-line arguments"""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "bamfile",
        type=pysam.AlignmentFile,
        help="sam/bam from which to calculate coverage",
    )
    return parser.parse_args()


def main():
    """main method of program"""
    args = parse_args()

    for ref_name, ref_length in zip(args.bamfile.references, args.bamfile.lengths):
        total_aligned_length = 0
        for segment in args.bamfile.fetch(ref_name):
            if not (
                segment.is_secondary or segment.is_supplementary or segment.is_unmapped
            ):
                total_aligned_length += segment.query_alignment_length
        print(
            "\t".join(
                [ref_name, str(total_aligned_length / ref_length), str(ref_length),]
            )
        )


if __name__ == "__main__":
    main()
