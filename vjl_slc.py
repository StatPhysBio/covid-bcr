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

def ham_dist_mat(strings):
    dists = np.zeros((len(strings), len(strings)))
    for i,s1 in enumerate(strings):
        for j,s2 in enumerate(strings):
            if i < j:
                break
            dists[i][j] = dists[j][i] = (hamming_distance(s1,s2)/len(s1))
    return dists

#  For making histogram to determine SLC threshold
def min_hams(mat):
    mat_copy = np.copy(mat)
    np.fill_diagonal(discopy, np.inf)
    return np.amin(discopy,axis=0)

def slc_dist_input(dist_mat, threshold):
    #  Convert distance matrix to 1d condensed distance matrix
    cdm = squareform(dist_mat)
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
    dists = ham_dist_mat(unique_cdr3s)

    _, clusters = slc_dist_input(dists, threshold)

    #  Map sequences to clusters
    clone_dict = {}
    for i, c in enumerate(clusters):
        clone_dict.setdefault(c, []).extend(seq_map[unique_cdr3s[i]])

    #  Save clutsers to lineage bin
    return clone_dict

def slc_all_bins(vjl_dict, partis=None, abstar=None, threshold=None):
    print(partis)
    lineage_dict = {}
    if not threshold:
        threshold = 0.15
    for i,key in enumerate(vjl_dict):
        lineage_dict[key] = slc_vjl_bin(vjl_dict[key], partis=partis, abstar=abstar, threshold=threshold)
    return lineage_dict

def vjl_slc(anns, partis=None, abstar=None, threshold=None):
    lineages = slc_all_bins(bin_vjl(anns, partis=partis, abstar=abstar), partis=partis, abstar=abstar, threshold=threshold)
    return lineages
