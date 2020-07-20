#!/usr/bin/env python3
"""
Trim N's from the beginnings and ends of all contigs
"""

import argparse
import textwrap
import re
import sys

import screed

starting_Ns = re.compile("^N+")
ending_Ns = re.compile("N+$")


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


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "fasta", type=screed.open, help="fasta file " "containing sequences to trim"
    )
    return parser.parse_args()


def trim_Ns(sequence_to_trim):
    """
    Remove all N's from beginning and end of a string

    Args:
        sequence_to_trim (str): sequence to remove N's from

    Returns:
        trimmed_seq (str): sequence_to_trim with N's removed
        num_trimmed_bases (int): number of N's removed
    """
    trimmed_seq = starting_Ns.sub("", sequence_to_trim)
    trimmed_seq = ending_Ns.sub("", trimmed_seq)
    return trimmed_seq, len(sequence_to_trim) - len(trimmed_seq)


def main():
    args = parse_args()

    for seq in args.fasta:
        if seq.sequence.startswith("N") or seq.sequence.endswith("N"):
            trimmed_seq, num_trimmed_bases = trim_Ns(seq.sequence)
            print(format_fasta(seq.name, trimmed_seq))
            print("\t".join([seq.name, str(num_trimmed_bases)]), file=sys.stderr)
        else:
            print(format_fasta(seq.name, seq.sequence))


if __name__ == "__main__":
    main()
