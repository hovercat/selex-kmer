#!/usr/bin/env python3

import argparse
import numpy as np

args_parser = argparse.ArgumentParser()
args_parser.add_argument("-i", type=str, required=True)
args_parser.add_argument("--cleaned", type=str, required=False)
args_parser.add_argument("-k", type=int, default=6)
args_parser.add_argument("-u", "--unique", action="store_true", default=False,
                         help="Unique-parameter: Count every kmer once per aptamer.")

import logging

alphabet_dna = {
    'A': 0,
    'C': 1,
    'G': 2,
    'T': 3
}

alphabet_dna_i = {
    0: 'A',
    1: 'C',
    2: 'G',
    3: 'T'
}


def sequence_to_int(seq):
    seq_hash = 0
    for i in range(0, len(seq)):
        seq_hash = seq_hash << 0b10  # one nucleotide is exactly two bits.
        seq_hash = seq_hash | alphabet_dna[seq[i]]

    return seq_hash


def int_to_seq(seq_hash, k):
    seq = []
    for i in range(0, k):
        seq.append(alphabet_dna_i[seq_hash & 0b11])
        seq_hash = seq_hash >> 2  # shift by 2 out
    return ''.join(reversed(seq))


def count_kmers(fasta_file_name, k, asv_fasta, uniques=False):
    # init kmer counter dict with 0 to 4 to pow k
    kmer_map = dict()
    for i in range(0, pow(4, k)):
        kmer_map[i] = 0

    # put ASVs into asv_list
    asv_list = list()
    with open(asv_fasta, 'r') as asv_fasta_file:
        for id in asv_fasta_file:
            seq = asv_fasta_file.readline().rstrip('\n')
            #counts = np.array(id.rstrip().split(' ')[1].split('-'), dtype=int)
            #total_count = np.sum(counts)
            #asv_list[seq] = counts[-1]  # hahahaha
            asv_list.append(seq)

    seq_list = list()
    with open(fasta_file_name, 'r') as fasta_file:
        for seq_id in fasta_file:
            seq = fasta_file.readline().rstrip('\n')

            # count kmers only if seq is not yet seen
            # only count kmers if seq in ASVs
            if seq in seq_list or seq not in asv_list:
                continue

            # ASV added to seq_dict
            seq_list.append(seq)

            kmer_list = [sequence_to_int(seq[i:i+k]) for i in range(0, len(seq) - k)]

            # if uniques is set, we do not count duplicate kmers in an aptamer e.g. kmer AAA in AAAA would appear two times.
            if uniques:
                kmer_list = list(set(kmer_list))  # removes duplicates


            for kmer_id in kmer_list:
                kmer_map[kmer_id] += 1

    return kmer_map


def print_csv(kmer_map, k):
    print("kmer;count\n", end="")
    for i in range(0, pow(4, k)):
        print("{};{}\n".format(int_to_seq(i, k), kmer_map[i]), end="")


def main():
    args = args_parser.parse_args()

    kmer_counts = count_kmers(args.i, args.k, args.cleaned, args.unique)
    print_csv(kmer_counts, args.k)


main()
