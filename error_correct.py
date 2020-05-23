from Bio.SeqIO import parse
import numpy as np
from typing import List, Tuple, TextIO
from jellyfish import hamming_distance

def get_cprimer(uid: str) -> str:
    return uid.split("|")[1].split("=")[-1]

def get_vprimer(uid: str) -> str:
    return uid.split("|")[2].split("=")[-1]

def get_abundance(uid: str) -> int:
    return int(uid.split("|")[3].split("=")[-1])

def get_time(uid: str) -> int:
    return int(uid.split("|")[4].split("=")[-1])

def get_replicate(uid: str) -> int:
    return int(uid.split("|")[6].split("=")[-1])

def fasta_parse(fasta: str) -> Tuple[List[str], List[str]]:
    uids = []
    seqs = []
    for seq in parse(fasta, 'fasta'):
        seqs.append(str(seq.seq))
        uids.append(seq.id)
    return uids, seqs

def write_to_fasta(save_name: str, headers: List[str], sequences: List[str]) -> None:
    with open(save_name, "w") as new_fasta:
        for i, header in enumerate(headers):
            new_fasta.write(">"+header + "\n" + sequences[i] + "\n")

def update_uid(uid: str, abundance: int) -> str:
    uid_split = uid.split("|")
    uid_split[3] = "DUPCOUNT=" + str(abundance)
    updated_uid = "|".join(uid_split)
    return updated_uid

def sort_data_by_abundance(uids: List[str], seqs: List[str]) -> Tuple[List[str], List[str]]:
    sorted_uids = sorted(uids, key=lambda u: get_abundance(u), reverse=True)
    sorted_seqs = [seq
                   for _,seq in sorted(zip(uids,seqs),
                                       key=lambda pair: get_abundance(pair[0]),
                                       reverse=True)]
    return sorted_uids, sorted_seqs

def error_correct_marginal(uids: List[str], seqs: List[str],
                           delta_r: np.float32 = 1.0,
                           delta_a: np.float32 = 1.0) -> Tuple[List[str], List[str]]:
    sorted_uids, sorted_seqs = sort_data_by_abundance(uids, seqs)
    sorted_abundances = [get_abundance(u) for u in sorted_uids]
    parent_indices = np.where(np.array(sorted_abundances) >= 10 ** (1 / delta_r))[0]

    for i in parent_indices:
        if sorted_abundances[i] == 0:
            continue
        for j,seq2 in reversed(list(enumerate(sorted_seqs))):
            if i == j:
                break
            if (sorted_abundances[i] == sorted_abundances[j]
                or sorted_abundances[j] == 0):
                continue
            abundance_log_ratio = (np.log10(sorted_abundances[i])
                                   - np.log10(sorted_abundances[j]))
            if 1 / abundance_log_ratio > delta_r:
                break
            if sorted_abundances[j] / abundance_log_ratio > delta_a:
                break
            seq1 = sorted_seqs[i]
            d = hamming_distance(seq1, seq2)
            if d / abundance_log_ratio <= delta_r:
                sorted_abundances[i] += sorted_abundances[j]
                sorted_abundances[j] = 0

    #  Get remaining sequences and updated uids
    remaining_seqs = []
    updated_uids = []

    for index, abund in enumerate(sorted_abundances):
        if abund == 0:
            continue
        remaining_seqs.append(sorted_seqs[index])
        updated_uids.append(update_uid(sorted_uids[index], sorted_abundances[index]))
    return updated_uids, remaining_seqs

def error_correct_total(uids: List[str], seqs: List[str],
                  d_tol: int = 1, a_tol: np.float32 = None) -> Tuple[List[str], List[str]]:
    if a_tol is None:
        a_tol = 1.0

    #  Sort by descending abundance
    sorted_uids, sorted_seqs = sort_data_by_abundance(uids, seqs)
    sorted_uids = sorted(uids, key=lambda u: get_abundance(u), reverse=True)
    sorted_seqs = [seq
                   for _,seq in sorted(zip(uids,seqs),
                                       key=lambda pair: get_abundance(pair[0]),
                                       reverse=True)]
    sorted_abundances = [get_abundance(u) for u in sorted_uids]

    for i,seq1 in enumerate(sorted_seqs):
        if sorted_abundances[i] == 0:
            continue
        for j,seq2 in reversed(list(enumerate(sorted_seqs))):
            if i == j:
                break
            if sorted_abundances[j] == 0:
                continue
            abundance_ratio = sorted_abundances[i] / sorted_abundances[j]
            if abundance_ratio < a_tol:
                break
            if hamming_distance(seq1, seq2) <= d_tol:
                sorted_abundances[i] += sorted_abundances[j]
                sorted_abundances[j] = 0

    #  Get remaining sequences and updated uids
    remaining_seqs = []
    updated_uids = []

    for index, abund in enumerate(sorted_abundances):
        if abund == 0:
            continue
        remaining_seqs.append(sorted_seqs[index])
        updated_uids.append(update_uid(sorted_uids[index], sorted_abundances[index]))
    return updated_uids, remaining_seqs

def group_data(uids: List[str], seqs: List[str]):
    #  Group by sequence length, c primer, v primer,
    #  time, and replicate
    num_header_info = len(uids[0].split("|"))
    if num_header_info == 7:
        key_func=lambda h_in,s_in: (len(s_in), get_cprimer(h_in),
                                    get_vprimer(h_in), get_time(h_in),
                                    get_replicate(h_in))
    elif num_header_info >= 5:
        key_func=lambda h_in,s_in: (len(s_in), get_cprimer(h_in),
                                    get_vprimer(h_in), get_time(h_in))
    elif num_header_info == 4:
        key_func=lambda h_in,s_in: (len(s_in), get_cprimer(h_in),
                                    get_vprimer(h_in))
    elif num_header_info < 4:
        print("Missing information in header.\n",
              "Cannot group sequences correctly." )
        return

    grouped_data = {}
    for u,s in zip(uids, seqs):
        key = key_func(u,s)
        if key not in grouped_data:
            grouped_data[key] = {'uids': [], 'sequences': []}
        grouped_data[key]['uids'].append(u)
        grouped_data[key]['sequences'].append(s)

    return grouped_data

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
    grouped_data = group_data(uids, seqs)

    for this_pass in range(args.passes):
        for key in grouped_data:
            out_uids, out_seqs = error_correct(grouped_data[key]['uids'],
                                               grouped_data[key]['sequences'],
                                               args.d_tol, args.a_tol)
            #  Update dictionary so things don't have to be resorted
            grouped_data[key]['uids'] = out_uids
            grouped_data[key]['sequences'] = out_seqs

    #  Combine all uids and sequences
    out_uids = []
    out_seqs = []
    for key in grouped_data:
        out_uids += grouped_data[key]['uids']
        out_seqs += grouped_data[key]['sequences']

    write_to_fasta(args.fasta.replace("_collapsed-unique.fasta","_corrected.fasta"), out_uids, out_seqs)

if __name__ == '__main__':
    main()
