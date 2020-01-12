#!/usr/bin/env python3
"""
Given a list of gene symbols and an ensembl gtf annotation of a genome,
outputs the coordinates for each gene in the list.
"""

import argparse
import collections
import re

attribute_regex = re.compile('(.+) "(.+)"')

Feature = collections.namedtuple(
    'Feature',
    ['chrom', 'feature_type', 'start_position', 'end_position', 'attributes']
)


def genes_list_type(file_location):
    """
    Parses a file containing a list of gene symbols, one per line, and
    returns a list of the symbols as strings
    """
    gene_symbols = []
    with open(file_location, 'r') as genes_list_file:
        for line in genes_list_file:
            gene_symbols.append(line.strip())
    return gene_symbols

def parse_attributes(attributes_string):
    """
    Parse the attributes field in a gtf entry into a dictionary mapping
    attribute name to attribute value.

    >>> parse_attributes('gene_id "GG034.2"; gene_name "ESR1";')
    {'gene_id': 'GG034.2', 'gene_name': 'ESR1'}
    """
    attributes_dict = {}
    # skip the last split because they all end with ';'
    for attribute in attributes_string.split(';')[:-1]:
        match = attribute_regex.match(attribute.strip())
        attributes_dict[match.group(1)] = match.group(2)
    return attributes_dict


def parse_gtf_line(gtf_line):
    """
    Parse a string containing a single line of a gtf file into a
    Feature object
    """
    fields = gtf_line.strip().split('\t')
    chrom, feature_type = fields[0], fields[2]
    start_position, end_position = int(fields[3]), int(fields[4])
    attributes = parse_attributes(fields[8])
    return Feature(
        chrom,
        feature_type,
        start_position,
        end_position,
        attributes,
    )


def gtf_type(file_location):
    """
    Given the path to an ensembl-style gtf, parse it and yield each
    line as a Feature object
    """
    with open(file_location, 'r') as gtf_file:
        for line in gtf_file:
            if not line.startswith('#'):
                yield parse_gtf_line(line)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        'genes_list', type=genes_list_type,
        help='file containing a list of gene symbols, one per line')
    parser.add_argument(
        'gtf', type=gtf_type,
        help='path to ensembl-style gtf file containing coordinates of genes')
    return parser.parse_args()


def main():
    args = parse_args()

    # make a dictionary mapping gene symbol string to a Feature object
    # corresponding to that gene. We fill in all the keys now in order
    # to keep them in order, but set the values to None
    symbols_to_features = collections.OrderedDict(
        [(symbol, None) for symbol in args.genes_list]
    )

    # now, loop through the gtf entries, filling in the
    # symbols_to_features dictionary as we encounter the symbols
    for feature in args.gtf:
        if (feature.feature_type == 'gene'
                and 'gene_name' in feature.attributes
                and feature.attributes['gene_name'] in symbols_to_features):
            symbols_to_features[feature.attributes['gene_name']] = feature

    # finally, loop through the list of input symbols, printing out the
    # coordinates for each symbol
    for symbol, feature in symbols_to_features.items():
        if feature is None:  # symbol not found in gtf file
            print('{}\tNA'.format(symbol))
        else:
            print('{}\t{}:{}-{}'.format(
                symbol,
                feature.chrom,
                feature.start_position,
                feature.end_position,
            ))


if __name__ == '__main__':
    main()

