#!/usr/bin/env python3

import argparse
import gzip
import collections

"""
A group of samples. The 'name' field is some identifier
string for the group; the 'members' field is a list of
the name of each sample in the group.
"""
Group = collections.namedtuple('Group', ['name', 'members'])

def make_header_string(groups):
    """
    Given a list of Group objects, make a header string for
    the output that lists all the default field names
    followed by each group name.
    """
    default_fields = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER']
    group_name_fields = list(map(lambda g: g.name, groups))
    return '\t'.join(default_fields + group_name_fields)

def gt_index_from_format(format_string):
    """
    Given a format string, return the index of the GT
    field. E.g., if FORMAT='GT:GQ:DP:HQ', this function
    will return 0.
    """
    for i, field in enumerate(format_string.split(':')):
        if field == 'GT':
            return i
    # TODO raise error if none of fields is 'GT'

def calculate_mafs(groups, gt_index, individual_fields):
    """
    Calculate the minor allele frequency within each group
    given the list of groups, the index of the genotype
    field of the format strings, and a dict mapping
    individual name string to individual format string.

    Arguments:
    - groups: an ordered list of Group objects
    - gt_index: the index of the genotype field of the
      format string. E.g., if FORMAT='GT:GQ:DP:HQ', this
      will be = 0
    - individual_fields: an ordered list of the raw
      contents of each individual's field, e.g.,
      '0|1:1:51:51'

    Returns: a list of minor allele frequencies with each
             group, in the same order as in the list of
             groups.
    """
    mafs = []
    for group in groups:
        total_allele_count = 0
        minor_allele_count = 0
        for member in group.members:
            gt_field = individual_fields[member].split(':')[gt_index]
            genotypes = tuple(map(lambda i: int(i),
                gt_field.replace('/', '|').split('|')))
            total_allele_count += len(genotypes)
            minor_allele_count += sum(genotypes)
        mafs.append(minor_allele_count / total_allele_count)
    return mafs

def parse_groups_file(groups_file):
    """
    Parses a file listing groups of samples in
    tab-separated format. First field is name of group;
    subsequent fields are the names of individuals within
    that group. Returns a list of Group objects.
    """
    groups = []
    for line in groups_file:
        splits = line.split('\t')
        groups.append(Group(splits[0], splits[1:]))
    return groups

def parse_args():
    parser = argparse.ArgumentParser(description='Group samples in a vcf into '
            'groups and then calculate the minor allele frequency for each '
            'variant in each group.')
    parser.add_argument('vcf', help='The vcf file to analyze. Can be gzipped.')
    parser.add_argument('groups_file', type=argparse.FileType('r'),
            help='A tab-separated file listing groups, one per line. The '
            'first field on each line is the name of the group and all '
            'subsequent fields are the names of individuals in that group.')
    return parser.parse_args()

def main():
    args = parse_args()

    groups = parse_groups_file(args.groups_file)

    if args.vcf.endswith('.gz'):
        vcf_file = gzip.open(filename, 'r')
    else:
        vcf_file = open(filename, 'r')

    for line in vcf_file:
        if line.startswith('##'):
            # this is one of the many header lines we don't
            # care about
            pass
        elif line.startswith('#'):
            # TODO protect against crashing on files without this line
            # this is the one header line we do care about
            individuals = line.split('\t')[9:]
            print(make_header(individuals))
        else: # this is an actual entry
            splits = line.split('\t')
            gt_index = gt_index_from_format(splits[8])
            mafs = calculate_mafs(groups, gt_index, splits[9:])
            print('\t'.join(splits[:7] + list(map(mafs, lambda i: str(i)))))

if __name__ == '__main__':
    main()
