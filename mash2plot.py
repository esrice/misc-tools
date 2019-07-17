#!/usr/bin/env python3

"""
mash2plot.py -- take mashmap output filtered to only one reference chromosome
    and one query chromosome and plot the alignment.
"""

import sys
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

def main():
    with open(sys.argv[1], 'r') as f:
        for line in f:
            splits = line.strip().split()
            x1, x2 = list(map(lambda i: float(i)/1e6, splits[7:9]))
            y1, y2 = list(map(lambda i: float(i)/1e6, splits[2:4]))
            if splits[4] == '+':
                plt.plot([x1, x2], [y1, y2], color='blue', marker='o')
            else:
                plt.plot([x1, x2], [y2, y1], color='red', marker='o')

    plt.xlabel(splits[5])
    plt.ylabel(splits[0])
    plt.savefig('plot.pdf')

if __name__ == '__main__':
    main()
