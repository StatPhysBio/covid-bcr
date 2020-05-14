import numpy as np
import pickle
from jellyfish import hamming_distance
from Bio.SeqIO import parse
from typing import List, Tuple, TextIO
import progressbar
from scipy.spatial.distance import squareform

def get_cprimer(uid: str) -> str:
    return uid.split("|")[1].split("=")[-1]

def get_vprimer(uid: str) -> str:
    return uid.split("|")[2].split("=")[-1]

def get_abundance(uid: str) -> int:
    return int(uid.split("|")[3].split("=")[-1])

def get_time(uid: str) -> int:
    return int(uid.split("|")[4].split("=")[-1])

def fasta_parse(fasta: str, patient, timepoint, severity, singletons=True):
    seqs = []
    headers = []
    for seq in parse(fasta, 'fasta'):
        uid, cprimer, vprimer, abundance = (x for x in seq.id.split('|'))

        abundance_int = int(abundance.split("=")[-1])
        if not singletons:
            if abundance_int == 1:
                continue

        header_info = [uid+"-"+patient, cprimer, vprimer,abundance,
                       "TIME=" + str(timepoint), "SEVERITY=" + severity]
        header = "|".join(header_info)
        header = ">" + header
        headers.append(header)
        seqs.append(str(seq.seq))
    return seqs, headers

def ham_dist_oneform(strings):
    num_seqs = len(strings)
    num_entries = int((num_seqs**2 - num_seqs) / 2)
    dists = np.zeros(num_entries,dtype=np.uint16)
    index = 0
    for i,s1 in enumerate(strings):
        for j,s2 in enumerate(strings):
            if i <= j:
                break
            dists[index] = hamming_distance(s1,s2)
            index+=1
    return dists

#  For making histogram to determine SLC threshold
def min_hams(oneform):
    mat = squareform(oneform)
    np.fill_diagonal(mat, 9999)
    return np.amin(mat,axis=0)

def group_data(uids: List[str], seqs: List[str]):
    #  Group by sequence length, c primer, and v primer
    grouped_data = {}
    for u,s in zip(uids, seqs):
        cprimer = get_cprimer(u)
        vprimer = get_vprimer(u)
        time = get_time(u)
        slen = len(s)
        key = (slen, cprimer, vprimer, time)
        if key not in grouped_data:
            grouped_data[key] = {'uids': [], 'sequences': []}
        grouped_data[key]['uids'].append(u)
        grouped_data[key]['sequences'].append(s)

    return grouped_data

def pickle_save(pickle_file, contents):
    import pickle
    with open(pickle_file, 'wb') as handle:
        pickle.dump(contents, handle, protocol=pickle.HIGHEST_PROTOCOL)

def main():
    import argparse
    parser = argparse.ArgumentParser(
        description='FASTA hamming distance to nearest')
    parser.add_argument('--fasta', type=str, help='path to FASTA')
    args = parser.parse_args()

    seqs,headers = fasta_parse(args.fasta,'0','0','0')
    grouped = group_data(headers, seqs)
    nearestdist = []

    bar = progressbar.ProgressBar(max_value=len(grouped))
    for i,key in enumerate(grouped):
        oneform=ham_dist_oneform(grouped[key]['sequences'])
        mhams = min_hams(oneform)
        nearestdist += mhams.tolist()
        bar.update(i)
    bar.finish()
    pickle_save(args.fasta.replace(".fasta","_min_hams.pickle"),nearestdist)

if __name__=='__main__':
    main()
