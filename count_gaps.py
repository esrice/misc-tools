#!/usr/bin/env python3

import argparse
import re
from Bio import SeqIO


def parse_args():
    parser = argparse.ArgumentParser(
        description="Given a fasta containing "
        "a genome assembly, count the gaps in each sequence and output a "
        "tab-separated file in which each line is a sequence and the "
        "fields are sequence name; number of gaps; comma-separated list "
        "of lengths of ungapped subsequences."
    )
    parser.add_argument(
        "assembly_fasta",
        type=lambda f: SeqIO.parse(f, "fasta"),
        help="a fasta file containing a genome assembly to count gaps in",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    for record in args.assembly_fasta:
        contig_lengths = list(map(len, re.split("[nN]+", str(record.seq))))
        print(
            "\t".join(
                [
                    record.id,
                    str(len(contig_lengths) - 1),
                    ",".join(map(str, contig_lengths)),
                ]
            )
        )


if __name__ == "__main__":
    main()
