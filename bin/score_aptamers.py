#!/usr/bin/env python3

import argparse
import math

import pandas as pd
import numpy as np
import scipy.stats as st

args_parser = argparse.ArgumentParser()
args_parser.add_argument("-i", type=str, required=True)
args_parser.add_argument("-f", "--fasta", type=str, required=True)
args_parser.add_argument("-k", type=int, required=True)
args_parser.add_argument("-n1", type=str, required=True)
args_parser.add_argument("-n2", type=str, required=True)


def main():
    args = args_parser.parse_args()

    fisher_csv = pd.read_csv(args.i, delimiter=';')
    scoring_vector = fisher_csv.set_index("kmer")
    scoring_vector = scoring_vector["opposite_z_value"]
    scoring_vector = scoring_vector.to_dict()

    print("sequence;name1;name2;scores;max_shifting_5;min_shifting_5\n", end="")

    aptamers = list()
    with open(args.fasta, 'r') as aptamers_fasta:
        for id in aptamers_fasta:
            total_count = np.sum(np.array(id.rstrip().split(' ')[1].split('-'), dtype=int))
            seq = aptamers_fasta.readline().rstrip()

            #            if total_count < 100:
            #                continue

            aptamers.append(seq)

    for seq in aptamers:
        scores = dict()
        kmers = len(seq) - args.k
        for i in range(0, kmers):
            kmer = seq[i:i + args.k]

            if scoring_vector[kmer] < -20:
                scores[kmer] = -20
            elif scoring_vector[kmer] > 20:
                scores[kmer] = 20
            else:
                scores[kmer] = scoring_vector[kmer]

        scores = np.array(list(scores.values()))
        score = np.mean(scores)

        max_shifting = 0
        min_shifting = 0
        shifting_frame = 5
        for i in range(0, kmers-shifting_frame):
            shifting_score = np.mean(scores[i:i+shifting_frame])
            if shifting_score > max_shifting:
                max_shifting = shifting_score
            if shifting_score < min_shifting:
                min_shifting = shifting_score

        print(seq, end=";")
        print(args.n1, end=";")
        print(args.n2, end=";")
        print(score, end=";")
        print(max_shifting, end=";")
        print(min_shifting, end="\n")
    # print_csv(dict1, dict2, fisher_dict, z_dict, args.k)


main()
