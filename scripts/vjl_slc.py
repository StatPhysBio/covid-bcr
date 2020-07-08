#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Script to construct lineages from annotated sequences.
    Copyright (C) 2020 Montague, Zachary

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

from jellyfish import hamming_distance
import numpy as np
from scipy.cluster.hierarchy import linkage,fcluster,dendrogram
from scipy.spatial.distance import squareform

def bin_vjl(annotations: list, partis: bool = False, abstar: bool = True) -> dict:
    """Bins annotations by V gene, J gene, and CDR3 length.

    Parameters
    ----------
    annotations : list
        List of all productive or unproductive annotations.
    partis :  bool, optional
        Specifies if the annotations came from partis.
    abstar : bool, optional
        Specifies if the annotations came from abstar.

    Returns
    -------
    vjl_dict : dict
        Nested dictionary which contains bins of annotations grouped
        by V gene, J gene, and CDR3 length.
    """

    vjl_dict = {}
    if partis:
        for annotation in annotations:
            v = annotation["v_gene"].split("*")[0]
            j = annotation["j_gene"].split("*")[0]
            l = annotation["cdr3_length"]
            if v not in vjl_dict:
                vjl_dict[v] = {}
            if j not in vjl_dict[v]:
                vjl_dict[v][j] = {}
            if l not in vjl_dict[v][j]:
                vjl_dict[v][j][l] = []

            vjl_dict[v][j][l].append(annotation)
        return vjl_dict
    elif abstar:
        for annotation in annotations:
            v = annotation["v_gene"]["gene"]
            j = annotation["j_gene"]["gene"]
            l = len(annotation["junc_nt"])
            if v not in vjl_dict:
                vjl_dict[v] = {}
            if j not in vjl_dict[v]:
                vjl_dict[v][j] = {}
            if l not in vjl_dict[v][j]:
                vjl_dict[v][j][l] = []

            vjl_dict[v][j][l].append(annotation)
        return vjl_dict

    else:
        print("Need an option for the annotation input!")

def ham_dist_vectorform(strings: list) -> np.array:
    """Constructs the Hamming distance vector-form for a VJL grouping used for clustering.

    This function takes in a set of strings, the observed CDR3s in a VJL grouping,
    and computes upper triangle of the Hamming (normalized by length) square matrix.
    This vector-form is what is used as input for the single-linkage clustering
    algorithm given by scipy. Additionally, using this in tandem with np.float16
    decreases the memory usage substantially.

    Parameters
    ----------
    strings : list of strings

    Returns
    -------
    dists : np.array
        Vector-form of normalized Hamming distances given by np.float16 precision.
    """

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

def get_minimum_distances(vectorform: np.array) -> np.array:
    """Returns the minimum pairwise Hamming distances, useful for determing a clustering threshold.

    Converts the vector-form distances (upper-triangle of the Hamming distance matrix)
    to a square matrix and obtains the minimums in each column.

    Parameters
    ----------
    vectorform : np.array
        Array of values which form the upper triangle of the Hamming distance matrix.

    Returns
    -------
    np.array
        np.array of the minimum of each column in the Hamming distance matrix.
    """

    mat_copy = np.copy(squareform(vectorform))
    np.fill_diagonal(mat_copy, np.inf)
    return np.amin(mat_copy,axis=0)

def slc_dist_input(vectorform: np.array, threshold: float) -> (np.array, np.array):
    """Performs single-linkage clustering on the Hamming distance vector-form array.

    Parameters
    ----------
    vectorform : np.array
        Array of values which form the upper triangle of the Hamming distance matrix.
    threshold : float
        Distance threshold for the single-linkage clustering algorithm.

    Returns
    -------
    links : np.ndarray
        The hierarchical clustering encoded as a linkage matrix. Though not used
        in this analysis, this is included because it is used to illustrate
        the clusterings using dendrograms.
    cluster : np.ndarray
        An array which specifies what string goes into which cluster.
    """

    # Make links according to the minimum distances.
    links=linkage(vectorform, method='single')
    #  Cluster using the given threshold.
    clusters = fcluster(links, threshold, criterion='distance')
    return links, clusters

def slc_vjl_bin(vjl_bin: list, partis: bool = False, abstar: bool = True,
                threshold: float = None) -> dict:
    """Obtains the set of unique CDR3s from a VJL grouping and performs single-linkage clustering.

    Parameters
    ----------
    vjl_bin : list
        List of annotations which have the same V gene, J gene, and CDR3 length.
    partis :  bool, optional
        Specifies if the annotations came from partis.
    abstar : bool, optional
        Specifies if the annotations came from abstar.
    threshold : float, optional
        Specifies the distance threshold to be use for single-linkage clustering.

    Returns
    -------
    clone_dict : dict
        Dictionary of annotations which represents clustering annotations in a VJL grouping.
    """

    #  Set default threshold.
    if not threshold:
        threshold = 0.15

    #  Create map from unique CDR3s to annotations in vjl bin.
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

    #  If there is only one unique CDR3, no clustering is necessary,
    if len(unique_cdr3s) < 2:
        return {1: annotations_map[unique_cdr3s[0]]}

    # Create distance matrix.
    dists = ham_dist_vectorform(unique_cdr3s)

    #  Clusters is a list of integers which maps
    #  each unique CDR3 to a cluster.
    _, clusters = slc_dist_input(dists, threshold)

    #  Place the corresponding annotation into each cluster.
    clone_dict = {}
    for i, c in enumerate(clusters):
        clone_dict.setdefault(int(c), []).extend(annotations_map[unique_cdr3s[i]])

    return clone_dict

def slc_all_bins(vjl_dict: dict, partis: bool = False, abstar: bool = True,
                 threshold: float = None) -> dict:
    """Perform single-linkage clustering over all VJL bins.

    Parameters
    ----------
    vjl_bin : list
        List of annotations which have the same V gene, J gene, and CDR3 length.
    partis :  bool, optional
        Specifies if the annotations came from partis.
    abstar : bool, optional
        Specifies if the annotations came from abstar.
    threshold : float, optional
        Specifies the distance threshold to be use for single-linkage clustering.

    Returns
    -------
    lineages : dict
        Nested dictionary which contains cluters of annotations grouped
        by V gene, J gene, CDR3 lenght, and cluster ID.
    """

    lineages = {}
    if not threshold:
        threshold = 0.15
    for v in vjl_dict:
        if v not in lineages:
            lineages[v] = {}
        for j in vjl_dict[v]:
            if j not in lineages[v]:
                lineages[v][j] = {}
            for l in vjl_dict[v][j]:
                if l not in lineages[v][j]:
                    lineages[v][j][l] = {}
                lineages[v][j][l] = slc_vjl_bin(vjl_dict[v][j][l], partis=partis,
                                                abstar=abstar, threshold=threshold)
    return lineages

def vjl_slc(annotations, partis: bool = False, abstar: bool = True,
            threshold: float = None) -> dict:
    """Perform VJL-grounping and single-linkage clustering over all annotations.

    Parameters
    ----------
    annotations : list
        List of all productive or unproductive annotations.
    partis :  bool, optional
        Specifies if the annotations came from partis.
    abstar : bool, optional
        Specifies if the annotations came from abstar.
    threshold : float, optional
        Specifies the distance threshold to be use for single-linkage clustering.

    Returns
    -------
    lineages : dict
        Nested dictionary which contains cluters of annotations grouped
        by V gene, J gene, CDR3 lenght, and cluster ID.
    """

    lineages = slc_all_bins(bin_vjl(annotations, partis=partis, abstar=abstar),
                            partis=partis, abstar=abstar, threshold=threshold)
    return lineages
