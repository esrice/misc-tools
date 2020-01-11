#!/usr/bin/env python3
"""
Given a list of gene symbols and an ensembl gtf annotation of a genome,
outputs the coordinates for each gene in the list.
"""

import argparse
import collections
import os.path
import re

attribute_regex = re.compile('(.+) "(.+)"')

def genes_list_type(file_location):
    """
    Parses a file containing a list of gene symbols, one per line, and
    returns a list of the symbols as strings
    """
    gene_symbols = []
    with open(file_location, 'r') as genes_list_file:
        for line in genes_list_file:
            gene_symbols.append(line.strip())


def parse_attributes(attributes_string):
    attributes_dict = {}
    # skip the last split because they all end with ';'
    for attribute in attributes_string.split(';')[:-1]:
        match = attribute_regex.match(attribute)
        attributes_dict[match.group(1)] =  match.group(2)
    return attributes_dict


def parse_gtf_line(gtf_line):
    fields = line.strip().split('\t')
    chrom, feature_type = fields[0], fields[2]
    start_position, end_position = int(fields[3]), int(fields[4])
    attributes = parse_attributes(fields[8])


def gtf_type(file_location):
    """
    Given the path to an ensembl-style gtf, parse it and yield each
    line as a Feature object
    """
    with open(file_location, 'r') as gtf_file:
        for line in gtf_file:
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
    # loop through all the gene symbols
    for gene_symbol in args.genes_list:
        # look up each gene symbol in the gff and then print a
        # tab-separated line containing the gene symbol, chromosome,
        # and start/end coordinates
        gene = args.gff[gene_symbol]
        print('\t'.join(map(
            str, [gene_symbol, gene.seqid, gene.start, gene.end]
        )))


if __name__ == '__main__':
    main()

