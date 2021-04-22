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
            sequence_nodes[seq] = { Sequence_Node(seq, id, total_counts) }

    seq_x = filtered_results['seq_x'].to_list()
    seq_y = filtered_results['seq_y'].to_list()
    for i in range(len(seq_x)):
        seq_x_i = seq_x[i]
        seq_y_i = seq_y[i]

        sequence_nodes[seq_x_i].update(sequence_nodes[seq_y_i])
        sequence_nodes[seq_y_i] = sequence_nodes[seq_x_i]
        #sequence_nodes[row.seq_x].add_friend(sequence_nodes[row.seq_y])
        #sequence_nodes[row.seq_y].add_friend(sequence_nodes[row.seq_x])

    friends_sets = list()
    for friends_set in sequence_nodes.values():
        friends_sets.append(frozenset(friends_set))
    friends_sets = set(friends_sets)


   # seqs_to_remove = set()
    best_friends = dict()
    for friends_set in friends_sets:
        best_friend = None
        max_count = 0

        for friend in friends_set:
            if friend.count > max_count:
                best_friend = friend
                max_count = friend.count

        best_friends[best_friend.seq] = best_friend
       # seqs_to_remove.update(friends_set)

    for seq, seq_node in best_friends.items():
        print(seq_node.id)
        print(seq_node.seq)


main()
