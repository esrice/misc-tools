#!/usr/bin/env python3

import argparse


def parse_args():
    """Parse command-line arguments"""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "busco", type=argparse.FileType("r"), help="filtered BUSCO output table"
    )
    parser.add_argument("n", type=int, help="total number of BUSCOs")
    return parser.parse_args()


def main():
    """main method of program"""
    args = parse_args()
    # this dictionary has BUSCO id as key, and as values:
    #   if the BUSCO is "Complete", the string "Complete"
    #   if the BUSCO is "Duplicated", the number of copies present
    buscos = {}
    # read through the file and fill up the dictionary
    with open(args.busco, "r") as busco_file:
        for line in busco_file:
            if not line.starts_with("#"):
                splits = line.strip().split()
                busco, status = splits[:2]
                # if the BUSCO is complete, denote it as such
                if status == "Complete":
                    buscos[busco] = "Complete"
                # if it is marked as "Duplicated", keep track of how many
                # times we see it
                elif status == "Duplicated":
                    if busco not in buscos:
                        buscos[busco] = 0
                    buscos[busco] += 1

    num_complete_buscos, num_duplicated_buscos = 0, 0
    for value in buscos.values():
        # if a BUSCO is marked as "Duplicated" but only appears once, it's
        # because we removed one or more of the instances, so it's actually
        # a complete BUSCO.
        if value == "Complete" or len(value) == 1:
            num_complete_buscos += 1
        else:
            num_duplicated_buscos += 1

    print(
        "Complete: {} ({}%)".format(
            num_complete_buscos, num_complete_buscos / args.n * 100
        )
    )
    print(
        "Duplicated: {} ({}%)".format(
            num_duplicated_buscos, num_duplicated_buscos / args.n * 100
        )
    )


if __name__ == "__main__":
    main()
