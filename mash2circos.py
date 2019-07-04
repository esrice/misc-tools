#!/usr/bin/env python3

import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='')
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

def mash_to_karyotype(mash_file, karyotype_outfile, query_color='orange',
        reference_color='green'):
    # read through the mashmap output to build dicts of
    # sequence names to sequence lengths for both the
    # query assembly and the reference assembly
    query_seq_sizes = {}
    reference_seq_sizes = {}
    for line in mash_file:
        splits = line.split(' ')
        query_name = splits[0]
        query_len = int(splits[1])
        reference_name = splits[5]
        reference_len = int(splits[6])

        if not query_name in query_seq_sizes:
            query_seq_sizes[query_name] = query_len
        if not reference_name in reference_seq_sizes:
            reference_seq_sizes[reference_name] = reference_len

    for query_name in sorted(query_seq_sizes.keys(), key=chr_sort_key):
        query_length = query_seq_sizes[query_name]
        print(' '.join(map(str, ['chr', '-', query_name,
            query_name.lstrip('chr'), 1, query_length, query_color])),
            file=karyotype_outfile)
    for reference_name in sorted(reference_seq_sizes.keys(), key=chr_sort_key):
        reference_length = reference_seq_sizes[reference_name]
        print(' '.join(map(str, ['chr', '-', reference_name,
            reference_name.lstrip('chr'), 1, reference_length,
            reference_color])), file=karyotype_outfile)

def main():
    args = parse_args()
    karyotype_file = open(args.output_prefix + '.karyotype', 'w')
    mash_to_karyotype(args.mashmap, karyotype_file)

if __name__ == "__main__":
    main()
