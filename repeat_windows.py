#!/usr/bin/env python3

import argparse
from Bio import SeqIO


def parse_args():
    parser = argparse.ArgumentParser(
        description="Given a fasta containing "
        "a genome assembly, divide it into windows and count the number "
        "of bases soft-masked as repetitive in each window, then output a "
        "file in bed-like format with columns: chr, start, end, fraction "
        "of total sequence in this window that is repetitive."
    )
    parser.add_argument(
        "assembly_fasta",
        type=lambda f: SeqIO.parse(f, "fasta"),
        help="a fasta file containing a genome assembly to count gaps in",
    )
    parser.add_argument(
        "-w",
        "--window-size",
        type=int,
        default=50000,
        help="size of windows to calculate repeat percentages in",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    for record in args.assembly_fasta:
        for window_start in range(0, len(record.seq), args.window_size):
            window_seq = record.seq[window_start : window_start + args.window_size]
            repeat_base_count = sum([int(c.islower()) for c in window_seq])
            repeat_fraction = repeat_base_count / len(window_seq)
            print(
                "\t".join(
                    map(
                        str,
                        [
                            record.id,
                            window_start + 1,
                            window_start + len(window_seq) + 1,
                            repeat_fraction,
                        ],
                    )
                )
            )


if __name__ == "__main__":
    main()
