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


def main():
    args = parse_args()

    # make list of ANN sub-field names
    ann_field_names = list(map(str.strip, annotation_re.match(
        vcf.infos['ANN'].desc).group(1).split('|')))

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
        ann_fields = dict(zip(
            ann_field_names,
            record.INFO['ANN'][0].split('|')
        ))
        output_fields += [ann_fields[n] for n in [
            'Gene_Name',
            'Annotation',
            'gene_region',  # TODO
            'effect',  # TODO
            'transcript_name',  # TODO
            'exon_number',  # TODO
            'HGVS.c',
            'HGVS.p',
        ]]


if __name__ == '__main__':
    main()
