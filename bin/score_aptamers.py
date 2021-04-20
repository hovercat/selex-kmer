#!/usr/bin/env python3

import argparse
import pandas as pd
import numpy as np


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
    scoring_vector = scoring_vector["p_value"]
    scoring_vector = scoring_vector.to_dict()

    print("sequence;name1;name2;score\n", end="")
    seq_dict = dict()


    aptamers = list()
    with open(args.fasta, 'r') as aptamers_fasta:
        for id in aptamers_fasta:
            seq = aptamers_fasta.readline().rstrip()
            aptamers.append(seq)

    kmer_scores = np.array(list(scoring_vector.values()))
    #modif = np.mean(kmer_scores[kmer_scores < 0.05])

    for seq in aptamers:
        scores = dict()
        kmers = len(seq)-args.k
        for i in range(0, kmers):
            kmer = seq[i:i+args.k]
            if 'GGG' in kmer or kmer.endswith('GG') or kmer.startswith('GG'):
                continue

            #if scoring_vector[kmer] < 0.2:
            scores[kmer] = scoring_vector[kmer]
            #if scores[kmer] < pow(10, -20):
            #    scores[kmer] = pow(10, -20)

            #scores.append(scoring_vector[kmer])
           # if scoring_vector[kmer] < 0.005:
           #     sig += 1
           # else:
           #     sig = 0

           # if sig > max_sig:
           #     max_sig = sig

            #if kmer_score > 200: #float("inf"):
            #    kmer_score = 200
            #elif kmer_score < -200: # float("inf")*-1:
            #    kmer_score = -200

           # scores[kmer] = kmer_score #scoring_vector[seq[i:i+args.k]]
            #scores.append(kmer_score)
        scores = list(scores.values())
#        seq_dict[seq] = score / kmers
#        scores = list(set(scores))
#        scores = np.sort((list(scores.values())))
#        scores = np.sort(scores)
#        scores = scores / np.mean(scores)
     ######   seq_dict[seq] = np.mean(scores[-10:])
        #seq_dict[seq] = np.sum(np.sign(scores)) * np.var(scores)
        #seq_dict[seq] = np.mean(scores) * np.std(scores)
        seq_dict[seq] = np.sum(np.array(scores) <  0.05) # modif) # / fisher_aptamer.iloc[j]["p_value"]
         #seq_dict[seq] = np.mean(scores) * np.var(scores / np.mean(scores))
        #scores = scores / np.mean(scores)
        #seq_dict[seq] = np.median(scores[-10:])

        #seq_dict[seq] = np.sum(scores) #/ kmers

#        seq_dict[seq] = max_sig / np.var(np.divide(scores, np.mean(scores)))
        #seq_dict[seq] = np.sum(scores)

#        seq_dict[seq] = np.mean(scores) * np.var(scores)

#        seq_dict[seq] = np.var(scores)
#        seq_dict[seq] = np.mean(scores[-10:]) - np.mean(scores[:10]) # doesnt make sense
        

        print(seq, end=";")
        print(args.n1, end=";")
        print(args.n2, end=";")
        print(seq_dict[seq], end="\n")
    #print_csv(dict1, dict2, fisher_dict, z_dict, args.k)


main()
