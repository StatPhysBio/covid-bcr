#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Script to cluster single cell sequences into lineages.
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

from abstar_pipeline import *
from vjl_slc import *
from utils import *
from jellyfish import hamming_distance
import pandas as pd
import numpy as np

def make_v_ref_dict():
    """Create a dictionary of V genes which end before the C beginning the CDR3
       specified by the anchor files.

    Parameters
    ----------
    None

    Returns
    -------
    v_dict : dict
        Dictionary where the key is the v gene and the item is the V gene sequence.
    """

    v_anchors = (pd.read_csv('/gscratch/stf/zachmon/covid/covid-bcr/sonia_input/V_gene_CDR3_anchors.csv')
                 .assign(g = v_anchors['gene'].str.split("*").str[0]))
    anchor_idx = v_anchors.drop_duplicates('g').set_index('g').loc[v_anchors['g'].unique()]['anchor_index']
    v_seqs, headers = fasta_read('/gscratch/stf/zachmon/covid/covid-bcr/igor_input/abstar_genomic_Vs.fasta')
    v_dict = {}
    for idx, h in enumerate(headers):
        gene = h.split("*")[0]
        if gene in v_dict:
            continue
        v_trim = v_seqs[idx][:anchor_idx.loc[gene]]
        #  Some V gene anchors are erroneous and end on the CDR3 rather than right before?
        if v_trim[-3:] == 'TGT':
            v_trim = v_trim[:-3]
        v_dict[gene] = v_trim
    v_dict = sort_dict(v_dict)
    return v_dict

def get_similar_v_genes():
    """Returns a dictionary of V genes 90% similar to a given V gene.

    Parameters
    ----------
    None

    Returns
    -------
    v_to_include : dict
        Dictionary where the keys are V genes and the items are V genes
        90% similar to key.
    """

    v_ref = make_v_ref_dict()
    v_ref_genes = list(v_ref.keys())
    v_ham_mat = np.zeros(shape=(len(v_ref), len(v_ref)))

    for idx1, v1 in enumerate(v_ref):
        for idx2, v2 in enumerate(v_ref):
            seq1 = v_ref[v1]
            seq2 = v_ref[v2]
            min_len = np.min([len(seq1), len(seq2)])
            #  Go backwards from where the CDR3 begins.
            v_ham_mat[idx1, idx2] = hamming_distance(seq1[-min_len:], seq2[-min_len:]) / min_len

    v_to_include = {}
    for idx1, v1 in enumerate(v_ref_genes):
        v_to_include[v1] = []
        for idx2, v2 in enumerate(v_ref_genes):
            if v_ham_mat[idx1, idx2] <= 0.1:
                v_to_include[v1].append(v2)
    return v_to_include

#  Make the dictionary a global variable so it's not recalculated unnecessarily.
v_to_include = get_similar_v_genes()

def cluster_singlecell_nt(singlecell_annotation, lineages, thresh):
    """Performs single-linkage clustering on a single cell given 90% V gene similarity.

    Parameters
    ----------
    singlecell_annotation : dict
        Dictionary containing information about the annotated single cell sequence.
    lineages : dict
        Unnested dictionary of clustered annotations [(V, J, L, cluster_id)].
    threshold : float
        Distance threshold for the single-linkage clustering algorithm.

    Returns
    -------
    to_insert : np.array
        np.array of keys of lineages into which single cells clustered successfully.
    """

    v_gene = singlecell_annotation['v_gene']['gene']
    cdr3 = singlecell_annotation['junc_nt']
    len_cdr3 = len(cdr3)

    subset = [key for key in lineages
              if key[2] == str(len_cdr3)
              and key[0] in v_to_include[v_gene]]

    min_distances = np.ones(len(subset))
    for idx, vjlc in enumerate(subset):
        distances = np.zeros(len(lineages[vjlc]), dtype=np.float16)
        for idxa, annotation in enumerate(lineages[vjlc]):
            distances[idxa] = hamming_distance(cdr3, annotation['junc_nt'])
        distances /= len_cdr3
        min_distances[idx] = (np.min(distances))

    to_insert = np.array(subset)[min_distances <= thresh]
    return to_insert

def cluster_singlecell_aa(singlecell_annotation, lineages, thresh):
    """Performs single-linkage clustering (using aa CDR3) on a single cell given 90% V gene similarity.

    Parameters
    ----------
    singlecell_annotation : dict
        Dictionary containing information about the annotated single cell sequence.
    lineages : dict
        Unnested dictionary of clustered annotations [(V, J, L, cluster_id)].
    threshold : float
        Distance threshold for the single-linkage clustering algorithm.

    Returns
    -------
    to_insert : np.array
        np.array of keys of lineages into which single cells clustered successfully.
    """

    v_gene = singlecell_annotation['v_gene']['gene']
    cdr3 = singlecell_annotation['junc_nt']
    cdr3_aa = translate(cdr3)
    len_cdr3 = len(cdr3)
    len_cdr3_aa = len(cdr3_aa)

    subset = [key for key in lineages
              if key[2] == str(len_cdr3)
              and key[0] in v_to_include[v_gene]]

    min_distances = np.ones(len(subset))
    for idx, vjlc in enumerate(subset):
        distances = np.zeros(len(lineages[vjlc]), dtype=np.float16)
        for idxa, annotation in enumerate(lineages[vjlc]):
            distances[idxa] = hamming_distance(cdr3_aa, translate(annotation['junc_nt']))
        distances /= len_cdr3_aa
        min_distances[idx] = (np.min(distances))

    to_insert = np.array(subset)[min_distances <= thresh]
    return to_insert

def main():
    import argparse
    parser = argparse.ArgumentParser(
        description='Clusters single cell sequences into lineages constructed previously '
                    'from bulk+plasma repertoire data using 90% v gene similarity and '
                    'single-linkage clustering. The output is saved as the input name but '
                    'replaces ".json" with "_with_sc.json".')
    parser.add_argument('--infile', type=str,
                        help='a .json lineage file.')
    parser.add_argument('--savedir', type=str,
                        help='path to where the output will be saved.')
    args = parser.parse_args()

    singlecell_annotations = json_open('/gscratch/stf/zachmon/covid/monoclonal_ntd_rbd.json')
    lineages = json_open(args.infile)
    denested_lineages = denest_lineages(lins['productive'])

    for annotation in singlecell_annotations:
        lins_clustered = cluster_singlecell_nt(annotation, denested_lineages, 0.15)
        for vjlc in lins_clustered:
            lineages['productive'][vjlc[0]][vjlc[1]][vjlc[2]][vjlc[3]] += [ann]

    save_name = args.infile.split('/')[-1].replace('.json','_with_sc.json')
    json_save(args.savedir + save_name, lineages)

if __name__ == '__main__':
    main()
