#!/usr/bin/env python3

import argparse
import math

args_parser = argparse.ArgumentParser()
args_parser.add_argument("-i1", type=str, required=True)
args_parser.add_argument("-i2", type=str, required=True)
args_parser.add_argument("-k", type=int, required=True)

import scipy.stats as stats
import csv
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

    with open(args.i1) as file1, open(args.i2) as file2:
        dict1 = kmer_csv_to_dict(file1)
        dict2 = kmer_csv_to_dict(file2)

    total1 = sum(dict1.values())
    total2 = sum(dict2.values())

    fisher_dict = dict()
    z_dict = dict()

    #print("kmer;count_x;count_y;p_value;odds;log2odds;z_value;opposite_z_value\n", end="")
    print("kmer;count_x;count_y;p_value;opposite_z_value\n", end="")
    for i in range(0, pow(4, args.k)):
        kmer = utils.int_to_seq(i, args.k)
        fisher_table = [[dict2.get(kmer, 0), total2], [dict1.get(kmer, 0), total1]]
        #fisher = stats.fisher_exact(fisher_table, alternative='less')
        fisher = stats.fisher_exact(fisher_table, alternative='greater')
        z = stats.norm.ppf(fisher[1])

        print(kmer, end=";")
        print(dict1.get(kmer, 0), end=";")
        print(dict2.get(kmer, 0), end=";")
        print(fisher[1], end=";")
        #print(fisher[0], end=";")
        #if fisher[0] > 0:
        #    print(math.log2(fisher[0]), end=";")
        #else:
        #    print("-inf", end=';')
        #print(z, end=";")
        print(z*(-1), end="\n")

    #print_csv(dict1, dict2, fisher_dict, z_dict, args.k)


main()
