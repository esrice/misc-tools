#!/usr/bin/env python3

import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='Make a circos plot from '
            'mashmap alignments.')
    parser.add_argument('mashmap', type=argparse.FileType('r'),
            help='output of mashmap')
    parser.add_argument('-p', '--output-prefix', default='out',
            help='prefix for output files [out]')
    return parser.parse_args()

# TODO for now, I'm just assuming that we're only making
# circos plots of sequences assigned to chromosomes with
# chr prefixes in their names because it makes life so
# much simpler. The next two functions are the beginnings
# of my feeble efforts to figure out how to generalize
# this, so I'm leaving them here in case I need them later
def sort_chromosome_names(chromosome_names):
    pass

def sort_sequence_names(sequence_names):
    """
    Sort names of sequences based on these rules:
    - all chromosomes should go before all unplaced scaffolds
    - chromosomes can be named with or without 'chr' prefix
    - sort chromosomes in numerical order followed by W, X, Y, Z, M
    - sort unplaced scaffolds in lexographical order
    """
    # this lambda returns true if argument is a chromosome name
    chromosome_lambda = lambda s: str.starts_with(s, 'chr') or s.isdigit() \
            or s in 'WXYZM'

    chrom_names = filter(chromosome_lambda, sequence_names)
    scaffold_names = filter(lambda s: not chromosome_lambda(x), sequence_names)

    return sort_chromosome_names(chrom_names) + sorted(scaffold_names)

def chr_sort_key(sequence_name):
    """
    Provides sort key for a chromosome such that chromosome
    names are sorted in order chr1, chr2, ..., chrN, chrX,
    chrY, chrM ('chr' prefixes optional)
    """
    extra_chroms_dict = {'W': 101, 'X': 102, 'Y': 103, 'Z': 104, 'M': 105}
    stripped_sequence_name = sequence_name.lstrip('chr')
    if stripped_sequence_name.isdigit():
        return int(stripped_sequence_name)
    elif stripped_sequence_name in extra_chroms_dict:
        return extra_chroms_dict[stripped_sequence_name]
    else:
        return 1000

def mash_to_circos(mash_file, karyotype_outfile, links_outfile,
        query_color='orange', reference_color='green'):
    """
    Given a mashmap output file, use it to make the
    karyotype and links input files to mashmap.

    Arguments:
    * mash_file: readable File pointing to mashmap output
    * karyotype_outfile: writable File to put lists of
      reference and query sequences in
    * links_outfile: writeable File to put alignment links
      from mashmap alignment into
    * query_color, reference_color: colors to use for
      differentiating query sequences from reference
      sequences in circos plot
    """
    # read through the mashmap output to build dicts of
    # sequence names to sequence lengths for both the
    # query assembly and the reference assembly
    query_seq_sizes = {}
    reference_seq_sizes = {}
    for line in mash_file:
        splits = line.strip().split(' ')
        query_name, query_len = splits[0], int(splits[1])
        query_start, query_end = int(splits[2]), int(splits[3])
        reference_name, reference_len = splits[5], int(splits[6])
        reference_start, reference_end = int(splits[7]), int(splits[8])

        if not query_name in query_seq_sizes:
            query_seq_sizes[query_name] = query_len
        if not reference_name in reference_seq_sizes:
            reference_seq_sizes[reference_name] = reference_len

        print(' '.join(map(str, [query_name + '_qry', query_start, query_end,
            reference_name + '_ref', reference_start, reference_end])),
            file=links_outfile)

    # print a line for each chromosome in the query and each in the ref
    for query_name in sorted(query_seq_sizes.keys(), key=chr_sort_key):
        query_length = query_seq_sizes[query_name]
        print(' '.join(map(str, ['chr', '-', query_name + '_qry',
            query_name.lstrip('chr'), 1, query_length, query_color])),
            file=karyotype_outfile)
    for reference_name in sorted(reference_seq_sizes.keys(), key=chr_sort_key,
            reverse=True):
        reference_length = reference_seq_sizes[reference_name]
        print(' '.join(map(str, ['chr', '-', reference_name + '_ref',
            reference_name.lstrip('chr'), 1, reference_length,
            reference_color])), file=karyotype_outfile)

def main():
    args = parse_args()
    karyotype_file = open(args.output_prefix + '.karyotype', 'w')
    links_file = open(args.output_prefix + '.links', 'w')
    mash_to_circos(args.mashmap, karyotype_file, links_file)

if __name__ == "__main__":
    main()
