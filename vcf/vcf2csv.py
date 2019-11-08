#!/usr/bin/env python3

import argparse
import re
import sys
import vcf

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

annotation_re = re.compile("Functional annotations: '(.+)'")

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('vcf', help='The vcf file to convert to csv',
                        type=lambda s: vcf.Reader(filename=s))
    return parser.parse_args()

def comma_join(fields):
    """
    Converts everything in the list to strings and then joins
    them with commas.
    """
    return ','.join(map(str,fields))

def make_annotation_dict(description):
    """
    Takes the 'desc' string of the ANN info header, as constructed by
    SNPeff, and builds a dictionary mapping each field name to its
    index.
    """
    annotation_dict = {}
    match = annotation_re.match(description)
    for i, field in map(str.strip, match.group(1).split('|')):
        annotation_dict[field] = i
    return annotation_dict

def main():
    args = parse_args()

    # make dictionary of ANN field name to index
    annotation_field_names = make_annotation_dict(vcf.infos['ANN'].desc)

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
        ann_fields = record.INFO['ANN'].split('|')
        output_fields += [
            ann_fields[ann_field_names['Gene_Name']],
            'seq_ontology',
            'gene_region',
            'effect',
            'transcript_name',
            'exon_number',
            'hgvs_c',
            'hgvs_p',
        ]


if __name__ == '__main__':
    main()
