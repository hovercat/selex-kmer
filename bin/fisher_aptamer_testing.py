#!/usr/bin/env python3

import argparse
import math

args_parser = argparse.ArgumentParser()
args_parser.add_argument("-i", type=str, required=True)
args_parser.add_argument("-r1", type=str, required=True)
args_parser.add_argument("-r2", type=str, required=True)
import scipy.stats as stats
import numpy as np
import utils
import pandas


def main():
    args = args_parser.parse_args()

    seqs = pandas.read_csv(args.i, delimiter='\t')
    seqs = seqs.set_index("sequence")
    totals = np.sum(seqs).to_dict()
    seqs = seqs.T.to_dict()


    print("seq;round_x;round_y;p_value\n", end="")

    r1 = args.r1
    r2 = args.r2

    for seq, seq_counts in seqs.items():
        if seq_counts[r1] == 0 and seq_counts[r2] == 0:
            p = 1.0
        else:
            p = stats.fisher_exact(
                [[seq_counts[r1],totals[r1]],
                 [seq_counts[r2],totals[r2]]],
                alternative="less"
            )[1]

        print(seq, end=";")
        print(r1, end=";")
        print(r2, end=";")
        print(p, end="\n")

main()
