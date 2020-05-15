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

#  Upper triangle of distance matrix
def ham_dist_oneform(strings):
    normalization = len(strings[0])
    num_seqs = len(strings)
    num_entries = int((num_seqs**2 - num_seqs) / 2)
    dists = np.zeros(num_entries,dtype=np.float16)
    index = 0
    for i,s1 in enumerate(strings):
        for s2 in strings[i+1:]:
            dists[index] = (hamming_distance(s1,s2) / normalization)
            index+=1
    return dists

#  For making histogram to determine SLC threshold
def min_hams(mat):
    mat_copy = np.copy(squareform(mat))
    np.fill_diagonal(mat_copy, np.inf)
    return np.amin(mat_copy,axis=0)

def slc_dist_input(dist_oneform, threshold):
    # Make links
    links=linkage(dist_oneform, method='single')
    #  Cluster according to threshold
    clusters = fcluster(links, threshold, criterion='distance')
    return links, clusters

def slc_vjl_bin(vjl_bin, partis=None, abstar=None, threshold=None):
    if not threshold:
        threshold = 0.15

    #  Create map from unique cdr3s to annotations in vjl bin
    if partis:
        annotations_map = {}
        for ann in vjl_bin:
            cdr3 = ann['cdr3_seqs'][0]
            if cdr3 not in annotations_map:
                annotations_map[cdr3] = []
            annotations_map[cdr3].append(ann)
    elif abstar:
        annotations_map = {}
        for ann in vjl_bin:
            cdr3 = ann['junc_nt']
            if cdr3 not in annotations_map:
                annotations_map[cdr3] = []
            annotations_map[cdr3].append(ann)
    else:
        print("Need an option for the annotation input!")

    unique_cdr3s = list(annotations_map.keys())
    if len(unique_cdr3s) < 2:
        return {1: annotations_map[unique_cdr3s[0]]}

    # Create distance matrix
    dists = ham_dist_oneform(unique_cdr3s)

    #  Clusters is a list of integers which maps
    #  each unique cdr3 to a cluster
    _, clusters = slc_dist_input(dists, threshold)

    #  Put annotations in clusters
    clone_dict = {}
    for i, c in enumerate(clusters):
        clone_dict.setdefault(c, []).extend(annotations_map[unique_cdr3s[i]])

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
