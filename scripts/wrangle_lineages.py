#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Script to wrangle data for analyses.
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

import numpy as np
import pandas as pd
import sys

from abstar_pipeline import (denest_lineages, merge_replicates_within_lineage,
                             get_lineage_progenitor_cdr3, translate)
from utils import *

def create_count_dict() -> dict:
    """Creates dictionary that is used to store counts.

    Parameters
    ----------
    None

    Returns
    -------
    count_dict : dict
        Dictionary used to store counts and single cell sequences.
    """

    orig_keys = ['unique', 'unique_merged', 'abundance', 'singleton']
    keys = ['bulk_' + k for k in orig_keys]
    keys += ['plasma_' + k for k in orig_keys]
    keys += ['combined_' + k for k in orig_keys
             if k != 'unique' and k != 'abundance']
    keys += ['rbd_unique', 'rbd_seqs', 'ntd_unique','ntd_seqs']
    count_dict = {}
    for k in keys:
        count_dict[k] = 0
        if 'seqs' in k:
            count_dict[k] = []
    return count_dict

def get_counts_by_time_replicate(headers: list, headers_unique_counts: list,
                                 headers_plasma_indistinct: list, times_and_replicates: list) -> dict:
    """Obtains unique, abundance, and singleton counts per replicate and time in a lineage.

    Parameters
    ----------
    headers : list
        List of sequence information which contains time and abundance counts.
    headers_unique_counts : list
        List of sequence information after replicate merging.
    headers_plasma_indistinct : list
        List of sequence information after replicate merging with plasma
        sequences can merge with bulk sequences.
    times_and_replicates : list
        List of combinations of times and replicates at which sequences in the
        lineage could have existed.

    Returns
    -------
    count_dict : dict
        Dictionary containing the unique, abundance, and singleton counts
        at all times and replicates in a lineage for the bulk and plasma
        repertoires in addition to single cell counts and sequences.
    lineage_primer : str
        The primer used by the majority of sequences in a lineage. (They should
        all have the same primer if the pipeline, which was supposed to remove
        sequences with mismatched primers, was run correctly.)
    """

    count_dict = create_count_dict()
    vprimers = []

    count_dict = {}
    for ti_rep in times_and_replicates:
        count_dict[ti_rep] = create_count_dict()

    for i,u in enumerate(headers_unique_counts):
        ti = get_time(u)
        rep = get_replicate(u)
        if 'plasma' in u:
            count_dict[(ti, rep)]['plasma_unique_merged'] += 1
        elif 'rbd' in u or 'ntd' in u:
            continue
        else:
            count_dict[(ti, rep)]['bulk_unique_merged'] += 1

    for i,u in enumerate(headers_plasma_indistinct):
        if 'rbd' in u or 'ntd' in u:
            continue
        vprimers.append(get_vprimer(u))
        ti = get_time(u)
        rep = get_replicate(u)
        abundance = get_abundance(u)
        count_dict[(ti,rep)]['combined_unique_merged'] += 1
        if abundance == 1:
            count_dict[(ti,rep)]['combined_singleton'] += 1
    lineage_primer = Counter(vprimers).most_common(1)[0][0]

    for u in headers:
        ti = get_time(u)
        rep = get_replicate(u)
        abundance = get_abundance(u)
        singleton = (abundance == 1)

        if 'plasma' in u:
            count_dict[(ti,rep)]['plasma_unique'] += 1
            count_dict[(ti,rep)]['plasma_abundance'] += abundance
            if singleton:
                count_dict[(ti,rep)]['plasma_singleton'] += 1
        elif 'rbd' in u:
            count_dict[(ti,rep)]['rbd_unique'] += 1
            count_dict[(ti,rep)]['rbd_seqs'].append(u.split('|')[0])
        elif 'ntd' in u:
            count_dict[(ti,rep)]['ntd_unique'] += 1
            count_dict[(ti,rep)]['ntd_seqs'].append(u.split('|')[0])
        else:
            count_dict[(ti,rep)]['bulk_unique'] += 1
            abundance = get_abundance(u)
            count_dict[(ti,rep)]['bulk_abundance'] += abundance
            if singleton:
                count_dict[(ti,rep)]['bulk_singleton'] += 1

    return count_dict,lineage_primer

def create_wrangled_df(in_lineages: dict, productive: bool = True) -> pd.DataFrame:
    """Obtains unique, abundance, and singleton counts per time (and time) for all lineages.

    Parameters
    ----------
    in_lineages : dict
        Nested dictionaries of lineages [V][J][L][cluster_id].
    productive : bool, optional
        Specify if wrangling productive lineages only.

    Returns
    -------
    df : pd.DataFrame
       DataFrame which specifies count and single cell information at each time
       and replicate for each lineage.
    """

    lineages = denest_lineages(in_lineages, productive=productive)
    patient = get_patient(lineages[list(lineages.keys())[0]][0]['seq_id'])

    times_and_replicates = []
    for vjlc in lineages:
        times_and_replicates += [(get_time(annotation['seq_id']),
                                  get_replicate(annotation['seq_id']))
                                 for annotation in lineages[vjlc]]
    times_and_replicates = sorted(list(set(times_and_replicates)),
                                  key = lambda k: (k[0], k[1]))
    num_ti_rep = len(times_and_replicates)

    dict_for_df = {'patient': [], 'v_gene': [], 'j_gene': [], 'cdr3_length': [], 'cluster_id': [],
                   'observed_cdr3_nt': [], 'observed_cdr3': [], 'progenitor_cdr3_nt': [],
                   'progenitor_cdr3': [], 'time': [], 'replicate':[], 'productive': [], 'primer': [],
                   'bulk_abundance': [], 'bulk_unique': [], 'bulk_unique_merged': [], 'bulk_singleton': [],
                   'plasma_abundance': [], 'plasma_unique': [], 'plasma_unique_merged': [], 'plasma_singleton': [],
                   'combined_unique_merged': [], 'combined_singleton': [],
                   'rbd_unique': [], 'rbd_seqs': [], 'ntd_unique': [], 'ntd_seqs': []}

    for vjlc in lineages:
        lineage = lineages[vjlc]
        lineage_unique_counts = merge_replicates_within_lineage(lineage)
        lineage_plasma_indistinct = merge_replicates_within_lineage(lineage, plasma_distinct=False)

        observed_cdr3_nt = Counter([ann['junc_nt'] for ann in lineage_unique_counts
                                    if 'monoclonal' not in ann['seq_id']]).most_common()[0][0]
        observed_cdr3_aa = translate(observed_cdr3_nt)
        progenitor_cdr3_nt = get_lineage_progenitor_cdr3(lineage_unique_counts, nt=True)
        progenitor_cdr3_aa = translate(progenitor_cdr3_nt)

        headers = [ann['seq_id'] for ann in lineage]
        headers_unique_counts = [ann['seq_id'] for ann in lineage_unique_counts]
        headers_plasma_indistinct = [ann['seq_id'] for ann in lineage_plasma_indistinct]

        count_dict,lineage_primer = get_counts_by_time_replicate(headers, headers_unique_counts,
                                                                 headers_plasma_indistinct, times_and_replicates)

        dict_for_df['patient'] += [patient] * num_ti_rep
        dict_for_df['v_gene'] += [vjlc[0]] * num_ti_rep
        dict_for_df['j_gene'] += [vjlc[1]] * num_ti_rep
        dict_for_df['cdr3_length'] += [vjlc[2]] * num_ti_rep
        dict_for_df['cluster_id'] += [vjlc[3]] * num_ti_rep
        dict_for_df['observed_cdr3_nt'] += [observed_cdr3_nt] * num_ti_rep
        dict_for_df['observed_cdr3'] += [observed_cdr3_aa] * num_ti_rep
        dict_for_df['progenitor_cdr3_nt'] += [progenitor_cdr3_nt] * num_ti_rep
        dict_for_df['progenitor_cdr3'] += [progenitor_cdr3_aa] * num_ti_rep
        dict_for_df['productive'] += [progenitor_cdr3_nt != ''] * num_ti_rep
        dict_for_df['primer'] += [lineage_primer] * num_ti_rep

        for ti_rep in count_dict:
            dict_for_df['time'].append(ti_rep[0])
            dict_for_df['replicate'].append(ti_rep[1])
            for count_type in count_dict[ti_rep]:
                dict_for_df[count_type].append(count_dict[ti_rep][count_type])

    df = pd.DataFrame(dict_for_df)

    return df

def main():
    import argparse

    parser = argparse.ArgumentParser(
        description='Wrangle lineage data for expansion, sharing, etc. analyses.')
    parser.add_argument('--lineages', type=str,
                        help='path to json lineages file')
    parser.add_argument('--productive', action='store_false', default=True,
                        help='wrangle observed productive lineages with '
                        'productive progenitors only (default: True)')
    parser.add_argument('--outfile', type=str,
                        help='path and name to save .csv file')
    args = parser.parse_args()

    if args.lineages is None:
        print('No input file')
        return
    if args.productive == True:
        print('Processing observed productive lineages with productive progenitors only.')
    else:
        print('Processing observed productive lineages with productive and unproductive progenitors.')

    lineages = json_open(args.lineages)

    #  Wrangle productive lineages.
    df_productive = create_wrangled_df(lineages['productive'], productive=args.productive)
    pd.DataFrame(df_productive).to_csv(args.outfile, index=False)
    print('Created', args.outfile)

    #  Wrangle unproductive lineages.
    df_unproductive = create_wrangled_df(lineages['unproductive'], productive=False)
    pd.DataFrame(df_unproductive).to_csv(args.outfile.replace('.csv', '_unproductive.csv'), index=False)
    print('Created',args.outfile.replace('.csv', '_unproductive.csv'))

if __name__ == '__main__':
    main()
