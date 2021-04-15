#!/usr/bin/env python3

import argparse

args_parser = argparse.ArgumentParser()
args_parser.add_argument("-i", type=str, required=True)
args_parser.add_argument("-k", type=int, default=6)
args_parser.add_argument("-u", "--unique", action="store_true", default=False, help="Unique-parameter: Count every kmer once per aptamer.")

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


def count_kmers(fasta_file_name, k, uniques=False):
    # init dict
    kmer_map = dict()
    for i in range(0, pow(4, k)):
        kmer_map[i] = 0
    
    with open(fasta_file_name, 'r') as fasta_file:
        for seq_id in fasta_file:
            seq = fasta_file.readline().rstrip('\n')
            kmer_list = list()
            for i in range(0, len(seq)-k):
                kmer_list.append(sequence_to_int(seq[i:i+k]))
            
            if uniques:
                kmer_list = list(set(kmer_list)) # removes duplicates
            
            for kmer_id in kmer_list:
                kmer_map[kmer_id] = kmer_map[kmer_id] + 1

    return kmer_map
    
def print_csv(kmer_map, k):
    print("kmer;count\n", end="")
    for i in range(0, pow(4,k)):
        print("{};{}\n".format(int_to_seq(i, k), kmer_map[i]), end="")

def main():
    args = args_parser.parse_args()

    kmer_counts = count_kmers(args.i, args.k, args.unique)
    print_csv(kmer_counts, args.k)
    
    
main()


