#!/usr/bin/env python3
"""
Given a vcf file containing deletions and a gff annotating the
reference genome, annotate each of the deletions as affecting either
an intergenic region, a regulatory region (+/- 2kb of a gene), an
intron, a coding sequence, or some combination of these.
"""

import argparse
import itertools
import os
import sys

import gffutils
import vcf


def gff_type(gff_path):
    """
    argparse type function for GFF files. Uses gffutils to create
    a database if one does not yet exist, and then returns that
    database.
    """
    db_path = gff_path + '.db'
    if gff_path.split('.')[-1] == 'db':
        return gffutils.FeatureDB(gff_path)
    elif os.path.exists(db_path):
        return gffutils.FeatureDB(db_path)
    else:
        print('Creating gff db...', file=sys.stderr)
        return gffutils.create_db(
            gff_path,
            db_path,
            # the id_spec is necessary because NCBI gff's do not follow
            # the GFF specification
            id_spec={'gene': 'db_xref'}
        )


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('vcf', help='the vcf file to annotate',
                        type=lambda f: vcf.Reader(filename=f))
    parser.add_argument('gff', help='annotations of the reference genome, '
                        'either in NCBI gff or a pre-created sqlite db',
                        type=gff_type)
    return parser.parse_args()


def get_deletion_effects(deletion_record, gff_db, regulatory_margin=2000):
    """
    Figure out the effects of a deletion using a gff.

    Args:
        deletion_record (vcf.Record): a vcf Record representing the
            deletion
        gff_db (gffutils.FeatureDB): a gffutils DB of a genome's
            annotations
        regulatory_margin (int): the amount of sequence on either side
            of the deletion to look in for genes to be classified as
            having their regulatory regions affected
    """
    affected_genes = set()
    intergenic = True
    regulatory, intronic, coding = [False] * 3

    # first, go through all the features that overlap the deletion
    # and use them to set the above booleans and add any affected
    # genes to affected_genes
    features_in_deletion = gff_db.region(
        seqid=deletion_record.CHROM,
        start=deletion_record.POS,
        end=deletion_record.sv_end
    )
    for feature in features_in_deletion:
        if feature.featuretype == 'gene':
            affected_genes.add(feature.attributes['Name'])
            intergenic = False
        elif feature.featuretype == 'intron':
            intronic = True
        elif feature.featuretype == 'CDS':
            coding = True

    # next, look for any genes *near* the deletion
    features_near_deletion = itertools.chain(
        gff_db.region(
            seqid=deletion_record.CHROM,
            start=deletion_record.POS - regulatory_margin,
            end=deletion_record.POS
        ),
        gff_db.region(
            seqid=deletion_record.CHROM,
            start=deletion_record.sv_end,
            end=deletion_record.sv_end + regulatory_margin,
        ),
    )
    for feature in features_near_deletion:
        if feature.featuretype == 'gene':
            regulatory = True
            intergenic = False
            affected_genes.add(feature.attributes['Name'])

    return affected_genes, intergenic, regulatory, intronic, coding


def annotate_deletion(record, affected_genes, intergenic, regulatory,
                      intronic, coding):
    """ adds INFO fields to a vcf record """
    # TODO fill this in
    return record


def main():
    """ __main__ method for this file """
    args = parse_args()

    writer = vcf.Writer(sys.stdout, args.vcf)
    # TODO fix output header to add new fields
    for record in args.vcf:
        effects_tuple = get_deletion_effects(record, args.gff,
                                             args.regulatory_margin)
        annotated_record = annotate_deletion(record, *effects_tuple)
        writer.write_record(annotated_record)
    writer.close()


if __name__ == '__main__':
    main()

