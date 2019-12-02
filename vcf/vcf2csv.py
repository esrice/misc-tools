#!/usr/bin/env python3

import argparse
import os
import re

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
    'effect',
    'hgvs_c',
    'hgvs_p',
]

annotation_re = re.compile("Functional annotations: '(.+)'")
default_ignore_annotations = [
    'intergenic_region',
    'intron_variant',
    '3_prime_UTR_variant',
    '5_prime_UTR_variant',
    'non_coding_transcript_exon_variant',
    'downstream_gene_variant',
    'upstream_gene_variant',
]


def comma_separated_list_type(cs_list):
    return cs_list.split(',')


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-o', '--output-dir', default='.',
                        help="directory where output should go ['.']")
    parser.add_argument(
        '-i', '--ignore-annotations', default=default_ignore_annotations,
        type=comma_separated_list_type,
        help='annotation types to exclude from output, separated by commas. '
        'Default: {}'.format(','.join(default_ignore_annotations)))
    parser.add_argument('vcf', help='The vcf file to convert to csv',
                        type=lambda s: vcf.Reader(filename=s))
    return parser.parse_args()


def comma_join(fields):
    """
    Converts everything in the list to strings and then joins
    them with commas.
    """
    return ','.join(map(str, fields))


def main():
    global output_field_names, annotation_re

    args = parse_args()

    # make list of ANN sub-field names
    ann_field_names = list(map(str.strip, annotation_re.match(
        args.vcf.infos['ANN'].desc).group(1).split('|')))

    # make one output csv file for each sample and write its header
    sample_outs = {}
    for sample_name in args.vcf.samples:
        sample_outs[sample_name] = open(os.path.join(
            args.output_dir, '{}.csv'.format(sample_name)), 'w')
        print(','.join(output_field_names + [sample_name]),
              file=sample_outs[sample_name])

    # make one joint output csv and write its header
    joint_outfile = open(os.path.join(args.output_dir, 'joint.csv'), 'w')
    print(','.join(output_field_names + args.vcf.samples),
          file=joint_outfile)

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
        ]))
        output_fields.append(record.INFO['AN'])
        ann_fields = dict(zip(
            ann_field_names,
            record.INFO['ANN'][0].split('|')
        ))

        # skip variants that are not interesting to the user
        if ann_fields['Annotation'] in args.ignore_annotations:
            continue

        output_fields += [ann_fields[n] for n in [
            'Gene_Name',
            'Annotation',
            'Annotation_Impact',
            'HGVS.c',
            'HGVS.p',
        ]]

        sample_fields = []
        for sample_name in args.vcf.samples:
            # write a line for any sample for which this variant is not 0/0
            genotype = record.genotype(sample_name)
            if genotype.is_variant:
                print(comma_join(output_fields + [genotype.data.GT]),
                      file=sample_outs[sample_name])

            # add this genotype to the list of sample genotypes to be
            # output jointly
            if genotype.called:
                sample_fields += genotype.gt_bases
            else:
                sample_fields += './.'

        # write a line for this variant in the joint output
        print(comma_join(output_fields + sample_fields), file=joint_outfile)


if __name__ == '__main__':
    main()
