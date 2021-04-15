#!/usr/bin/env python3

import argparse
import math
import pandas as pd

args_parser = argparse.ArgumentParser()
args_parser.add_argument("-i", type=str, required=True)
args_parser.add_argument("-f", "--fasta", type=str, required=True)
args_parser.add_argument("-k", type=int, required=True)

import scipy.stats as stats
import utils


def kmer_csv_to_dict(file, skip_header=True):
    if skip_header:
        file.readline()

    csv_dict = dict()
    for line in file:
        line = line.rstrip()
        key = line.split(';')[0]
        value = line.split(';')[1]
        csv_dict[key] = int(value)

    return csv_dict


def main():
    args = args_parser.parse_args()
    # args = args_parser.parse_args(["-i1", "output_EF07/kmers_R00.csv", "-i2", "output_EF07/kmers_R11.csv", "-k", 6])

    fisher_csv = pd.read_csv(args.i, delimiter=';')
    scoring_vector = fisher_csv.set_index("kmer")["log2odds"]
    scoring_vector = scoring_vector.to_dict()

    print("sequence;score\n", end="")
    seq_dict = dict()
    fasta = open(args.fasta)

    for id in fasta:
        id = id.rstrip()
        seq = fasta.readline().rstrip()

        score = 0
        kmers = len(seq)-args.k
        for i in range(0, kmers):
            score += scoring_vector[seq[i:i+args.k]]
        seq_dict[seq] = score / kmers

        print(seq, end=";")
        print(score, end="\n")
    #print_csv(dict1, dict2, fisher_dict, z_dict, args.k)


main()
