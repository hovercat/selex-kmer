
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
