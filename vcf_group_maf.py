#!/usr/bin/env python3

import argparse
import gzip
import collections
import yaml
import vcf
import sys

"""
A group of samples. The 'name' field is some identifier
string for the group; the 'members' field is a list of
the name of each sample in the group.
"""
Group = collections.namedtuple('Group', ['name', 'abbreviation', 'members'])

def parse_groups_yaml(groups_yaml):
    """
    Parses a file listing groups of samples in yaml format
    Returns a list of Group objects.
    """
    groups = []
    for group in yaml.load(groups_yaml)['groups']:
        try:
            name, abbreviation = group['name'], group['abbreviation']
            members = group['members']
        except:
            print('Groups yaml is missing a field. Aborting.',
                    file=sys.stderr)
            sys.exit(1)

        # TODO check to make sure abbreviation is not reserved
        groups.append(Group(name, abbreviation, members))

    return groups

def calc_group_maf(group, vcf_record):
    """
    Calculate the minor allele frequency for a given group
    of samples.

    Arguments:
    - group: a Group containing the samples to examine
    - vcf_record: a vcf Record containing the calls
    """
    # we're only dealing with biallelic sites for now
    if len(vcf_record.ALT) > 1:
        return '.'

    # the gt_type field is 0 for hom ref, 1 for het, 2 for
    # hom alt, and None for no call, so this is some super
    # convenient math
    calls = map(lambda s: vcf_record.genotype(s).gt_type, group.members)
    calls = list(filter(lambda c: c is not None, calls))

    if len(calls) == 0:
        return '.'
    else:
        return sum(calls)/(len(calls)*2)

def parse_args():
    parser = argparse.ArgumentParser(description='Group samples in a vcf into '
            'groups and then calculate the minor allele frequency for each '
            'variant in each group.')
    parser.add_argument('vcf', help='The vcf file to analyze. Can be gzipped.')
    parser.add_argument('groups_yaml', type=argparse.FileType('r'),
        help='A yaml file listing groups. See sample file for an example.')
    return parser.parse_args()

def main():
    args = parse_args()

    groups = parse_groups_yaml(args.groups_yaml)
    vcf_file = vcf.Reader(filename=args.vcf)

    # add INFO lines to header of vcf_file because it will
    # be used as a template for the output file
    for group in groups:
        vcf_file.infos[group.abbreviation] = vcf.parser._Info(
                group.abbreviation, 1, 'Float',
                'MAF of variant among group {}'.format(group.name), None, None)

    vcf_writer = vcf.Writer(sys.stdout, vcf_file)

    for record in vcf_file:
        if len(record.ALT) == 1: # only dealing with biallelic sites for now
            mafs = {g.abbreviation: calc_group_maf(g, record) for g in groups}
            record.INFO.update(mafs)
            vcf_writer.write_record(record)

if __name__ == '__main__':
    main()
