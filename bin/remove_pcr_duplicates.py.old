#!/usr/bin/env python3

import argparse

from Sequence_Node import Sequence_Node

args_parser = argparse.ArgumentParser()
args_parser.add_argument("-f", type=str, required=True)
args_parser.add_argument("-b", type=str, required=True)
args_parser.add_argument("--dereplicate", action="store_true", default=False)

import pandas as pd
import numpy as np


def main():
    args = args_parser.parse_args()

    blast = pd.read_csv(args.b,
                        delimiter="\t",
                        header=0,
                        names=[
                            "seq_x",
                            "seq_y",
                            "identity",
                            "length",
                            "mismatches",
                            "gaps",
                            "start_x",
                            "end_x",
                            "start_y",
                            "end_y",
                            "e_val",
                            "bitscore"
                        ])

    blast_filter = (blast.seq_x != blast.seq_y)
    filtered_results = blast[blast_filter]

    # read_sequence_file
    sequence_nodes = dict()
    with open(args.f, "r") as fasta_file:
        for id in fasta_file:
            id = id.rstrip()
            seq = fasta_file.readline().rstrip()
            total_counts = np.sum(np.array(id.rstrip().split(' ')[1].split('-'), dtype=int))
            sequence_nodes[seq] = Sequence_Node(seq, id, total_counts)

    for index, row in filtered_results.iterrows():
        sequence_nodes[row.seq_x].add_friend(sequence_nodes[row.seq_y])
        sequence_nodes[row.seq_y].add_friend(sequence_nodes[row.seq_x])

    seqs_to_remove = set()
    best_friends = dict()

    for seq, seq_node in sequence_nodes.items():
        if seq_node in seqs_to_remove or seq_node in best_friends:
            continue

        best_friend = seq_node.find_best_friend()
        best_friends[best_friend.seq] = best_friend
        removal_seqs = list(seq_node.friends)
        removal_seqs.remove(best_friend)
        seqs_to_remove.update(removal_seqs)

    for seq, seq_node in sequence_nodes.items():
        if seq in best_friends:
            print(seq_node.id)
            print(seq_node.seq)


main()
