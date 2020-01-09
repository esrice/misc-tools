#!/usr/bin/env python3
"""
Given a list of gene symbols and a gff annotation of a genome, outputs
the coordinates for each gene in the list.
"""

import argparse
import os.path

import gffutils


def genes_list_type(file_location):
    """
    Parses a file containing a list of gene symbols, one per line, and
    yields each gene symbol string
    """
    with open(file_location, 'r') as genes_list_file:
        for line in genes_list_file:
            yield line.strip()


def gff_type(file_location):
    """
    Given the path to a gff, check whether a gffutils database for this
    file exists yet. If so, open that database and return it; if not,
    create a new database at '{file_location}.symbol.db' indexed by
    gene symbol and return it.
    """
    db_path = '{}.symbol.db'.format(file_location)
    if os.path.exists(db_path):
        return gffutils.FeatureDB(db_path)
    else:
        # TODO this may require different options depending on which
        # reference versions we end up using, as nobody follows the
        # GFF specification correctly
        return gffutils.create_db(
            file_location,
            db_path,
            id_spec='symbol',
            keep_order=False,              # these two lines are not important
            sort_attribute_order=False,    # but speed up the database creation
        )


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        'genes_list', type=genes_list_type,
        help='file containing a list of gene symbols, one per line')
    parser.add_argument(
        'gff', type=gff_type,
        help='path to gff file containing coordinates of genes')
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

