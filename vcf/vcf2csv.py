#!/usr/bin/env python3

import argparse
import vcf
import sys

output_field_names = [
    'chromosome',
    'position',
    'ref',
    'alt',
    'alt_allele_frequency',
    'alt_allale_count',
    'num_alleles',
    'genes',
    'seq_ontology',
    'gene_region',
    'effect',
    'transcript_name',
    'exon_number',
    'hgvs_c',
    'hgvs_p',
]

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('vcf', help='The vcf file to convert to csv',
                        type=lambda s: vcf.Reader(filename=s))
    return parser.parse_args()

def comma_join(fields):
    return ','.join(fields)

def main():
    args = parse_args()

    # make dictionary of ANN field name to index
    vcf.infos['ANN'].desc

    output_field_names += args.vcf.samples
    print(','.join(output_field_names))  # header
    for record in args.vcf:
        output_fields = []
        output_fields += [
            record.CHROM,
            record.POS,
            record.REF,
        ]
        output_fields += list(map(comma_join, [
            record.ALT,
            record.INFO['AF'],
            record.INFO['AC'],
            record.INFO['AN'],
        ]))


if __name__ == '__main__':
    main()
