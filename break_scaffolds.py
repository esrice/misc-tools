#!/usr/bin/env python3

"""
break_scaffolds.py -- break scaffolds at gaps (i.e., strings of 'N')
"""

import argparse
import re
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO


def parse_args():
    parser = argparse.ArgumentParser(
        description="Given a fasta containing "
        "a genome assembly, break all sequences into contigs."
    )
    parser.add_argument(
        "assembly_fasta",
        type=lambda f: SeqIO.parse(f, "fasta"),
        help="a fasta file containing an assembly to break into contigs",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    for record in args.assembly_fasta:
        # TODO AFAIK, re.split can't return positions or
        # ranges or anything, so we'll need to do something
        # a little bit more complicated to add coordinates
        # to the sequence deflines
        for i, contig_seq in enumerate(re.split("[nN]+", str(record.seq))):
            contig_name = "{}.contig{}".format(record.id, i)
            contig_record = SeqRecord(Seq(contig_seq), id=contig_name, description="")
            print(contig_record.format("fasta"))


if __name__ == "__main__":
    main()
