#!/usr/bin/env python3
"""
Given a bed file containing the mapping coordinates of Hi-C reads, make
a heatmap.
"""

import argparse
import sys

import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-b', '--num-bins', type=int, default=500,
                        help='Number of bins to use for heatmap')
    parser.add_argument('-o', '--output', default='heat.png',
                        help='Place to put output plot [heat.png]')
    parser.add_argument('bed_file', type=argparse.FileType('r'),
                        help='Bed file containing the mapping coordinates of '
                        'Hi-C reads')
    parser.add_argument(
        'contigs_list', type=argparse.FileType('r'),
        help='list of contigs to make a heatmap of. Must be a tsv '
        'file with three columns: contig name; orientation; size in bp')
    return parser.parse_args()


def read_contigs_list(contigs_list_file):
    # contigs contains tuples with contig name, contig orientation, contig size,
    # and contig start coordinate
    contigs = []
    running_start = 0
    for line in contigs_list_file:
        name, orientation, size = line.strip().split()
        size = int(size)
        contigs.append((name, orientation, size, running_start))
        running_start += size

    return contigs


def contig_to_heat_coord(contig_name, contig_coord, contigs, contigs_index):
    """
    Convert from a contig:position coordinate to a single position.

    Arguments:
    contig_name: name of contig
    contig_coord: position on this contig
    contigs: DataFrame containing all the contigs
    contigs_index: dict mapping contig name to index in DataFrame
    """
    # figure out if this contig is in the index; if it isn't, return None.
    index = contigs_index.get(contig_name)
    if index is None:
        return None

    contig_row = contigs.loc[index]
    if contig_row['orientation'] == '+':
        return contig_row['start'] + contig_coord
    return contig_row['start'] + contig_row['size'] - contig_coord


def heat_coords_from_bed(bed_file, contigs, contigs_index):
    """
    Given a bed file and a list of contigs, make a list of coordinates
    representing read pair mapping locations.
    """
    # loop through the bed file
    heat_coords_x = []
    heat_coords_y = []
    last_contig, last_start, _, last_read_name = '', -1, -1, ''
    for line in bed_file:
        contig, start, end, read_name, _, _ = line.strip().split()
        start, end = int(start), int(end)

        # if the current read is a forward one, just store it
        # in anticipation of its reverse partner
        if read_name.endswith('/1'):
            last_contig, last_start = contig, start
            last_read_name = read_name
        elif read_name[:-2] == last_read_name[:-2]:
            last_coord = contig_to_heat_coord(last_contig, last_start,
                                              contigs, contigs_index)
            this_coord = contig_to_heat_coord(contig, start,
                                              contigs, contigs_index)

            if last_coord is not None and this_coord is not None:
                heat_coords_x.append(last_coord)
                heat_coords_y.append(this_coord)
                heat_coords_x.append(this_coord)
                heat_coords_y.append(last_coord)

    return heat_coords_x, heat_coords_y


def main():
    args = parse_args()

    # make a data frame with the contigs info in it
    contigs = pd.DataFrame(read_contigs_list(args.contigs_list),
            columns=['name', 'orientation', 'size', 'start'])

    # make a dictionary mapping the contig name to its index in the data frame
    contigs_index = { row['name']: i for i, row in contigs.iterrows() }

    heat_coords_x, heat_coords_y = heat_coords_from_bed(args.bed_file,
            contigs, contigs_index)

    fig, ax = plt.subplots()
    ax.hist2d(heat_coords_x, heat_coords_y, cmap=plt.cm.Reds,
            norm=mpl.colors.LogNorm(), bins=(args.num_bins, args.num_bins))
    ax.set_aspect('equal', 'box')
    for i, contig in contigs.iterrows():
        ax.axvline(x=contig['start'], linestyle='--', color='grey',
                linewidth=1)
        ax.axhline(y=contig['start'], linestyle='--', color='grey',
                linewidth=1)
    plt.xticks(contigs['start'], contigs['name'], rotation='vertical')
    plt.yticks(contigs['start'], contigs['name'])
    plt.savefig(args.output)

if __name__ == "__main__":
    main()
