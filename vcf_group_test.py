#!/usr/bin/env python3

import argparse
import gzip
import collections
import yaml
import vcf
import sys

def parse_args():
    parser = argparse.ArgumentParser(description='Filter vcf to only keep '
            'variants for which all or almost all of the individuals in a '
            'group of samples meets some requirement.')
    parser.add_argument('vcf', help='The vcf file to analyze. Can be gzipped.')
    parser.add_argument('samples', type=lambda s: s.split(','),
        help='a comma-separated list of samples')
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-a', '--hom-alt', action='store_true',
            help='members of group must be homozygous alternate for variant')
    group.add_argument('-r', '--hom-ref', action='store_true',
            help='memebrs of group must be homozygous reference for variant')
    group.add_argument('-z', '--het', action='store_true',
            help='members of group must be heterozygous for variant')
    parser.add_argument('-n', '--number-fails-allowed', type=int, default=0,
        help='number of individuals allowed to not comply with requirement')

    return parser.parse_args()

def main():
    args = parse_args()

    vcf_reader = vcf.Reader(filename=args.vcf)
    vcf_writer = vcf.Writer(sys.stdout, vcf_reader)

    if args.hom_ref:
        call_type = 0
    elif args.het:
        call_type = 1
    elif args.hom_alt:
        call_type = 2

    # only dealing with biallelic sites for now
    for record in filter(lambda r: len(r.ALT) == 1, vcf_reader):
        calls = map(lambda s: record.genotype(s).gt_type, args.samples)
        calls = filter(lambda c: c is not None, calls)
        pass_count = sum(map(lambda c: c == call_type, calls))
        if pass_count >= len(args.samples) - args.number_fails_allowed:
            vcf_writer.write_record(record)

if __name__ == '__main__':
    main()
