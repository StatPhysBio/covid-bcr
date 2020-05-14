from scipy.cluster.hierarchy import linkage,fcluster,dendrogram
from scipy.spatial.distance import squareform
import numpy as np
from jellyfish import hamming_distance

def bin_vjl(anns, partis=None, abstar=None):
    vjl_dict = {}
    if partis:
        for ann in anns:
            v = ann["v_gene"].split("*")[0]
            j = ann["j_gene"].split("*")[0]
            l = ann["cdr3_length"]
            key = (v,j,l)

            if key not in vjl_dict:
                vjl_dict[key] = []
            vjl_dict[key].append(ann)
        return vjl_dict
    elif abstar:
        for ann in anns:
            v = ann["v_gene"]["gene"]
            j = ann["j_gene"]["gene"]
            l = len(ann["junc_nt"])
            key = (v,j,l)

            if key not in vjl_dict:
                vjl_dict[key] = []
            vjl_dict[key].append(ann)
        return vjl_dict
    else:
        print("Need an option for the annotation input!")

def get_coarsegrained_bins(vjl_bins):
    v_bins = {}
    vj_bins = {}
    for vjl in vjl_bins:
        vj = (vjl[0],vjl[1])
        v = vjl[0]
        if v not in v_bins:
            v_bins[v] = []
        if vj not in vj_bins:
            vj_bins[vj] = []
        vj_bins[vj] += vjl_bins[vjl]
        v_bins[v] += vjl_bins[vjl]
    return v_bins, vj_bins

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
def min_hams(mat):
    mat_copy = np.copy(squareform(mat))
    np.fill_diagonal(mat_copy, np.inf)
    return np.amin(mat_copy,axis=0)

def slc_dist_input(dist_oneform, threshold):
    #  Convert distance matrix to 1d condensed distance matrix
    cdm = dist_oneform
    # Make links
    links=linkage(cdm, method='single')
    #  Cluster according to threshold
    clusters = fcluster(links, threshold, criterion='distance')

    return links, clusters

def slc_vjl_bin(vjl_bin, partis=None, abstar=None, threshold=None):
    if not threshold:
        threshold = 0.15

    #  Create map from unique cdr3s to sequences in vjl bin
    if partis:
        seq_map = {}
        for ann in vjl_bin:
            cdr3 = ann['cdr3_seqs'][0]
            if cdr3 not in seq_map:
                seq_map[cdr3] = []
            seq_map[cdr3].append(ann)
    elif abstar:
        seq_map = {}
        for ann in vjl_bin:
            cdr3 = ann['junc_nt']
            if cdr3 not in seq_map:
                seq_map[cdr3] = []
            seq_map[cdr3].append(ann)
    else:
        print("Need an option for the annotation input!")

    unique_cdr3s = list(seq_map.keys())
    if len(unique_cdr3s) < 2:
        return {1: seq_map[unique_cdr3s[0]]}

    # Create distance matrix
    dists = ham_dist_oneform(unique_cdr3s)

    _, clusters = slc_dist_input(dists, threshold)

    #  Map sequences to clusters
    clone_dict = {}
    for i, c in enumerate(clusters):
        clone_dict.setdefault(c, []).extend(seq_map[unique_cdr3s[i]])

    #  Save clusters to lineage bin
    return clone_dict

def slc_all_bins(vjl_dict, partis=None, abstar=None, threshold=None):
    print(partis)
    lineage_dict = {}
    if not threshold:
        threshold = 0.15
    print(threshold)
    for i,key in enumerate(vjl_dict):
        lineage_dict[key] = slc_vjl_bin(vjl_dict[key], partis=partis, abstar=abstar, threshold=threshold)
    return lineage_dict

def vjl_slc(anns, partis=None, abstar=None, threshold=None):
    lineages = slc_all_bins(bin_vjl(anns, partis=partis, abstar=abstar),
                            partis=partis, abstar=abstar, threshold=threshold)
    return lineages
