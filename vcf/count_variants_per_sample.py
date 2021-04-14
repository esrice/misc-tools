#!/usr/bin/env python3
"""
Count the number of called variants per sample in a VCF file.
"""

import argparse
import collections

import vcf


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-t",
        "--separate-by-type",
        help="separate counts by SV type",
        action="store_true",
    )
    parser.add_argument(
        "vcf", help="the vcf file to analyze", type=lambda f: vcf.Reader(filename=f)
    )
    return parser.parse_args()


def count_ignoring_sv_type(vcf_reader: vcf.Reader):
    """
    Count calls per sample, without taking SV type into account, and
    output a table.

    Args:
        vcf_reader(vcf.Reader): input vcf file
    """
    # these three Counters have sample name as key and number of
    # variants for that sample as values
    call_counts = collections.Counter()
    hom_alt_counts = collections.Counter()
    het_counts = collections.Counter()

    for record in filter(lambda r: not r.is_filtered, vcf_reader):
        for call in filter(lambda s: not s.is_filtered, record.samples):
            call_counts[call.sample] += 1
            if call.is_variant:
                if call.is_het:
                    het_counts[call.sample] += 1
                else:
                    hom_alt_counts[call.sample] += 1

    print("\t".join(["sample", "call_count", "hom_alt_count", "het_count"]))
    for sample in call_counts.keys():
        print(
            "\t".join(
                map(
                    str,
                    [
                        sample,
                        call_counts[sample],
                        hom_alt_counts[sample],
                        het_counts[sample],
                    ],
                )
            )
        )


def count_with_sv_type(vcf_reader: vcf.Reader):
    # these three defaultdicts have SV type as key, and a Counter as
    # value; the Counters have sample name as key, and variant
    # count as value. For example, call_counts['DEL']['sample1']
    # contains the number of deletions called for sample1.
    call_counts = collections.defaultdict(collections.Counter)
    hom_alt_counts = collections.defaultdict(collections.Counter)
    het_counts = collections.defaultdict(collections.Counter)

    for record in filter(lambda r: not r.is_filtered, vcf_reader):
        sv_type = record.INFO["SVTYPE"]
        for call in filter(lambda s: not s.is_filtered, record.samples):
            call_counts[sv_type][call.sample] += 1
            if call.is_variant:
                if call.is_het:
                    het_counts[sv_type][call.sample] += 1
                else:
                    hom_alt_counts[sv_type][call.sample] += 1

    print("\t".join(["sample", "sv_type", "call_count", "hom_alt_count", "het_count"]))
    for sv_type in call_counts.keys():
        for sample in call_counts[sv_type].keys():
            print(
                "\t".join(
                    map(
                        str,
                        [
                            sample,
                            sv_type,
                            call_counts[sv_type][sample],
                            hom_alt_counts[sv_type][sample],
                            het_counts[sv_type][sample],
                        ],
                    )
                )
            )


def main():
    args = parse_args()

    if args.separate_by_type:
        count_with_sv_type(args.vcf)
    else:
        count_ignoring_sv_type(args.vcf)


if __name__ == "__main__":
    main()
