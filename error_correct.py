from Bio.SeqIO import parse
import numpy as np
from typing import Tuple, TextIO
from jellyfish import hamming_distance

def get_abundance(uid):
    return int(uid.split("|")[3].split("=")[-1])

def fasta_parse(fasta: str):
    uids = []
    seqs = []
    for seq in parse(fasta, 'fasta'):
        seqs.append(str(seq.seq))
        uids.append(seq.id)
    return uids, seqs

def write_to_fasta(save_name: str, headers, sequences):
    with open(save_name, "w") as new_fasta:
        for i, header in enumerate(headers):
            new_fasta.write(header + "\n")
            new_fasta.write(sequences[i] + "\n")

def update_uid(uid: str, abundance: int):
    uid_split = uid.split("|")
    uid_split[3] = "DUPCOUNT=" + str(abundance)
    updated_uid = "|".join(uid_split)
    return updated_uid
def error_correct(uids, seqs, d_tol: int = 1, a_tol: np.float32 = None):
    if a_tol is None:
        a_tol = 1.0

    #  Sort by descending abundance
    sorted_uids = sorted(uids, key=lambda e: get_abundance(e), reverse=True)
    sorted_seqs = [seq
                   for _,seq in sorted(zip(uids,seqs),
                                       key=lambda pair: get_abundance(pair[0]),
                                       reverse=True)]
    sorted_abundances = [get_abundance(u) for u in sorted_uids]

    parent_child_dict = {}
    for i,seq1 in enumerate(sorted_seqs):
        if sorted_abundances[i] == 0:
            continue
        parent_child_dict[(i, seq1)] = []
        for j,seq2 in reversed(list(enumerate(sorted_seqs))):
            if i == j:
                continue

            if sorted_abundances[j] == 0:
                continue

            abundance_ratio = sorted_abundances[i] / sorted_abundances[j]
            if abundance_ratio < a_tol:
                break

            if hamming_distance(seq1, seq2) <= d_tol:
                parent_child_dict[(i, seq1)].append((j, seq2))
                sorted_abundances[i] += sorted_abundances[j]
                sorted_abundances[j] = 0

        if not parent_child_dict[(i, seq1)]:
            del parent_child_dict[(i, seq1)]

    #  Get remaining sequences and updated uids
    remaining_seqs = []
    updated_uids = []

    for index, abund in enumerate(sorted_abundances):
        if abund == 0:
            continue
        remaining_seqs.append(sorted_seqs[index])
        updated_uids.append(update_uid(sorted_uids[index], sorted_abundances[index]))
    return updated_uids, remaining_seqs

def main():
    """
    usage: python error_correct.py -h"""
    import argparse
    import pickle

    parser = argparse.ArgumentParser(
        description='FASTA error correction')
    parser.add_argument('fasta', type=str, help='path to FASTA')
    parser.add_argument('outbase', type=str, help='basename for output files')
    parser.add_argument('--d_tol', type=float, default=1,
                        help='Hamming distance tolerance (default 1)')
    parser.add_argument('--a_tol', type=float, default=None,
                        help='abundance ratio tolerance (if None, not '
                             'applied)')
    parser.add_argument('--passes', type=int, default=1,
                        help='number of times to repeat greedy clustering '
                             '(default 1)')
    parser.add_argument('--keep_singletons', action='store_true',
                        help="don't discard uncorrected singletons")
    args = parser.parse_args()

    uids, seqs = fasta_parse(args.fasta)
    for this_pass in range(1, args.passes + 1):
        uids, seqs = error_correct(uids, seqs, args.d_tol, args.a_tol)

    write_to_fasta(args.outbase + ".corrected.fa", uids, seqs)

if __name__ == '__main__':
    main()
