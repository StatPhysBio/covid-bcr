#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Script to get count information and run Fisher exact test.
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
import scipy.stats as stats

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
        Dictionary used to store counts.
    """

    orig_keys = ['unique', 'unique_merged', 'abundance', 'singleton']
    keys = orig_keys
    keys += ['plasma_' + k for k in orig_keys]
    keys += ['monoclonal_unique', 'monoclonal_seqs']
    count_dict = {}
    for k in keys:
        count_dict[k] = {}
    return count_dict

def get_counts_by_time_replicate(headers: list, headers_unique_counts: list, times: list, replicates: list) -> dict:
    """Obtains unique, abundance, and singleton counts per replicate and time in a lineage.

    Parameters
    ----------
    headers : list
        List of sequence information which contains time and abundance counts.
    headers_unique_counts : list
        List of sequence information after replicate merging.
    times : list
        List of times at which sequences in the lineage could have existed.

    Returns
    -------
    count_dict : dict
        Dictionary containing the unique, abundance, and singleton counts
        at all times and replicates in a lineage.
    lineage_primer : str
        The primer used by the majority of sequences in a lineage. (They should
        all have the same primer if the pipeline, which was supposed to remove
        sequences with mismatched primers, was run correctly.)
    """

    count_dict = create_count_dict()
    vprimers = [None]*len(headers_unique_counts)

    for ti in times:
        #  Instead of running a loop over all the data and obtaining
        #  the replicates, we will put 0's for all possible replicates.
        #  The replicates which don't exist will be removed quickly in get_counts.
        for r in replicates:
            for count_type in count_dict:
                count_dict[count_type][(ti,r)] = 0
                if count_type == 'monoclonal_seqs':
                    count_dict[count_type][(ti,r)] = []

    for i,u in enumerate(headers_unique_counts):
        vprimers[i] = get_vprimer(u)
        ti = get_time(u)
        rep = get_replicate(u)
        if 'plasma' in u:
            count_dict['plasma_unique_merged'][(ti,rep)] += 1
        elif 'monoclonal' in u:
            continue
        else:
            count_dict['unique_merged'][(ti,rep)] += 1

    lineage_primer = Counter(vprimers).most_common(1)[0][0]

    for u in headers:
        ti = get_time(u)
        rep = get_replicate(u)
        abundance = get_abundance(u)
        singleton = (abundance == 1)

        if 'plasma' in u:
            count_dict['plasma_unique'][(ti,rep)] += 1
            count_dict['plasma_abundance'][(ti,rep)] += abundance
            if singleton:
                count_dict['plasma_singleton'][(ti,rep)] += 1
        elif 'monoclonal' in u:
            count_dict['monoclonal_unique'][(ti,rep)] += 1
            count_dict['monoclonal_seqs'][(ti,rep)].append(u.split('|')[0])
        else:
            count_dict['unique'][(ti,rep)] += 1
            abundance = get_abundance(u)
            count_dict['abundance'][(ti,rep)] += abundance
            if singleton:
                count_dict['singleton'][(ti,rep)] += 1

    return count_dict,lineage_primer

def get_counts(in_lineages: dict, productive: bool = True,
               replicate: bool = False, abstar: bool = True, mergedcdr3=False) -> pd.DataFrame:
    """Obtains unique, abundance, and singleton counts per time (and time) for all lineages.

    Parameters
    ----------
    in_lineages : dict
        Nested dictionaries of lineages [V][J][L][cluster_id].
    patient_key : str
        Patient ID. This is used to acquire the times that will be present in a lineage
        without any unnecessary loops.
    replicate : bool, optional
        Split counts by replicate.
    abstar : bool, optional
        Specify if annotations came from abstar or partis.

    Returns
    -------
    df_counts : dict
        Dictionary of three pandas.Dataframes, one each for unique, abundance,
        and singletons counts.
    """

    primer_split_counts = create_count_dict()

    lineages = denest_lineages(in_lineages, productive=productive)
    patient_key = get_patient(lineages[list(lineages.keys())[0]][0]['seq_id'])

    replicates = []
    times = []
    for key in lineages:
        times += [get_time(ann['seq_id']) for ann in lineages[key]]
        replicates += [get_replicate(ann['seq_id']) for ann in lineages[key]]
    replicates = list(set(replicates))
    times = list(set(times))
    print(replicates)
    print(times)

    for key in lineages:
        lineage = lineages[key]
        lineage_unique_counts = merge_replicates_within_lineage(lineage)
        if not mergedcdr3:
            observed_cdr3_nt = Counter([ann['junc_nt'] for ann in lineage
                                        if 'monoclonal' not in ann['seq_id']]).most_common()[0][0]
            observed_cdr3_aa = translate(observed_cdr3_nt)
            progenitor_cdr3_nt = get_lineage_progenitor_cdr3(lineage, nt=True)
            progenitor_cdr3 = translate(progenitor_cdr3_nt)
        else:
            observed_cdr3_nt = Counter([ann['junc_nt'] for ann in lineage_unique_counts
                                        if 'monoclonal' not in ann['seq_id']]).most_common()[0][0]
            observed_cdr3_aa = translate(observed_cdr3_nt)
            progenitor_cdr3_nt = get_lineage_progenitor_cdr3(lineage_unique_counts, nt=True)
            progenitor_cdr3 = translate(progenitor_cdr3_nt)

        cdr3_list = [observed_cdr3_nt, observed_cdr3_aa,progenitor_cdr3_nt, progenitor_cdr3]
        dict_key = tuple(list(key) + cdr3_list)

        #  Annotations come from abstar.
        if abstar:
            headers = [ann['seq_id'] for ann in lineage]
            headers_unique_counts = [ann['seq_id'] for ann in lineage_unique_counts]
        #  Annotations come from partis.
        else:
            headers = [ann['unique_ids'][0] for ann in lineage]
            headers_unique_counts = [ann['unique_ids'][0] for ann in lineage_unique_counts]

        count_dict,lineage_primer = get_counts_by_time_replicate(headers, headers_unique_counts, times, replicates)
        if lineage_primer not in primer_split_counts['unique']:
            for count_type in primer_split_counts:
                primer_split_counts[count_type][lineage_primer] = {}
        for count_type in primer_split_counts:
            primer_split_counts[count_type][lineage_primer][dict_key] = count_dict[count_type]

    #  Convert to dataframe
    df_counts = {}
    for count_type in primer_split_counts:
        df_counts[count_type] = []
        #  Sort primer dictionary so primers are alphabetical.
        primer_split_counts[count_type] = sort_dict(primer_split_counts[count_type])
        for i, lineage_primer in enumerate(primer_split_counts[count_type]):
            df_counts[count_type].append(pd.DataFrame.from_dict(primer_split_counts[count_type][lineage_primer],
                                                                orient='index', dtype=np.uint32))

    for count_type in df_counts:
        #  Combine counts from all primers into a single dataframe
        df_counts[count_type] = pd.concat(df_counts[count_type], keys=primer_split_counts[count_type].keys())
        #  Remove columns of all zeros
        #df_counts[count_type] = df_counts[count_type].loc[:, (df_counts[count_type] != 0).any(axis=0)]
        #  Resetting indices causes the primer, v gene, j gene, cdr3 length, and cluster id
        #  to become columns. This is much more useful when using the fisher exact test
        #  for isolating the counts by primer.
        df_counts[count_type] = df_counts[count_type].reset_index()
        #  Rename the columns so that they are meaningful.
        df_counts[count_type] = df_counts[count_type].rename(columns={'level_0': 'primer',
                                                                      'level_1': 'v_gene',
                                                                      'level_2': 'j_gene',
                                                                      'level_3': 'cdr3_length',
                                                                      'level_4': 'cluster_id',
                                                                      'level_5': 'observed_cdr3_nt',
                                                                      'level_6': 'observed_cdr3',
                                                                      'level_7': 'progenitor_cdr3_nt',
                                                                      'level_8': 'progenitor_cdr3'})
    return df_counts

def get_columns_with_ints(df: pd.DataFrame) -> list:
    """Retrieves columns in the pd.DataFrame which have integers in them.

    Parameters
    ----------
    df : pd.DataFrame

    Returns
    -------
    cols_with_nums : list
        List of columns in df which are integers or are a tuple of integers.
    """
    cols_with_nums = []
    for col in df.columns:
        if type(col) == int:
            cols_with_nums.append(col)
        elif type(col[0]) == int:
            cols_with_nums.append(col)
    return cols_with_nums

def create_csv_for_analysis(lineages: dict, productive: bool = True,
                            replicate: bool = True, mergedcdr3: bool = False) -> dict:
    """Creates a dictionary from the count information and other information necessary for the R expansion analysis.

    Parameters
    ----------
    lineages : dict
        Nested dictionaries of clustered annotations [V][J][L][cluster_id].
    productive : bool, optional
        Specifies whether or not lineages empty-string CDR3s should be kept.
    replicate : bool, optional
        Specifies whether or not replicate informtion should be separated.

    Returns
    -------
    csv_dict : dict
        Dictionary of data to be used in R expansion analysis.
    """

    csv_dict = {'patient': [], 'v_gene': [], 'j_gene': [], 'cdr3_length': [],
                'observed_cdr3_nt': [], 'observed_cdr3': [], 'progenitor_cdr3_nt': [],
                'progenitor_cdr3': [], 'cluster_id':[], 'time': [], 'replicate':[],
                'abundance': [], 'unique': [], 'unique_merged': [], 'singleton': [],
                'productive': [],'primer':[], 'plasma_abundance': [], 'plasma_unique': [],
                'plasma_unique_merged': [], 'plasma_singleton': [], 'monoclonal_unique': [],
                'monoclonal_seqs': []}

    df_counts = get_counts(lineages, productive=productive,
                            abstar=True, mergedcdr3=mergedcdr3)
    ti_for_lins = get_columns_with_ints(df_counts['abundance'])
    lins = df_counts['abundance'][['v_gene','j_gene','cdr3_length','cluster_id','primer',
                                   'observed_cdr3_nt', 'observed_cdr3',
                                   'progenitor_cdr3_nt', 'progenitor_cdr3']].values.tolist()
    patient = None
    i = 0
    while patient is None:
        try:
            p = get_patient(lineages[lins[i][0]][lins[i][1]][lins[i][2]][lins[i][3]][0]['seq_id'])
            patient = int(p)
        except:
            i += 1
            pass

    for i, lin in enumerate(lins):
        for ti in ti_for_lins:
            csv_dict['primer'].append(lin[4])
            csv_dict['patient'].append(patient)
            csv_dict['v_gene'].append(lin[0])
            csv_dict['j_gene'].append(lin[1])
            csv_dict['cdr3_length'].append(int(lin[2])-6)
            csv_dict['observed_cdr3_nt'].append(lin[5])
            csv_dict['observed_cdr3'].append(lin[6])
            csv_dict['progenitor_cdr3_nt'].append(lin[7])
            csv_dict['progenitor_cdr3'].append(lin[8])
            csv_dict['cluster_id'].append(int(lin[3]))
            csv_dict['time'].append(ti[0])
            csv_dict['replicate'].append(ti[1])
            csv_dict['abundance'].append(df_counts['abundance'].iloc[i][ti])
            csv_dict['unique'].append(df_counts['unique'].iloc[i][ti])
            csv_dict['unique_merged'].append(df_counts['unique_merged'].iloc[i][ti])
            csv_dict['singleton'].append(df_counts['singleton'].iloc[i][ti])
            csv_dict['plasma_abundance'].append(df_counts['plasma_abundance'].iloc[i][ti])
            csv_dict['plasma_unique'].append(df_counts['plasma_unique'].iloc[i][ti])
            csv_dict['plasma_unique_merged'].append(df_counts['plasma_unique_merged'].iloc[i][ti])
            csv_dict['plasma_singleton'].append(df_counts['plasma_singleton'].iloc[i][ti])
            csv_dict['monoclonal_unique'].append(df_counts['monoclonal_unique'].iloc[i][ti])
            csv_dict['monoclonal_seqs'].append(df_counts['monoclonal_seqs'].iloc[i][ti])
            csv_dict['productive'].append(productive)
    return csv_dict

def make_csv_for_analysis(savename: str, csv_dict: dict) -> None:
    """Saves the necessary information to a csv for R expansion analysis.

    Parameters
    ----------
    savename : str
        Path and name of csv file to be saved.
    csv_dict : dict
        Dictionary of data necessary for expansion analysis.

    Returns
    -------
    None
    """

    col_list = ['patient','primer','productive',
                'v_gene','j_gene','cdr3_length','cluster_id',
                'common_cdr3','time','replicate',
                'abundance','unique','unique_merged','singleton',
                'plasma_abundance','plasma_unique','plasma_unique_merged','plasma_singleton',
                'monoclonal_unique', 'monoclonal_seqs']

    col_titles = ','.join(col_list)

    with open(savename,'w') as f:
        f.write(col_titles + '\n')
        for d_index,entry in enumerate(csv_dict['patient']):
            for index,key in enumerate(col_list):
                f.write(str(csv_dict[key][d_index]))
                if index != len(col_list) - 1:
                    f.write(',')
                else:
                    f.write('\n')

    print('Made', savename)

#  TODO: Finish Fisher exact test implementation here.
#        R implementation in notebook.

def fisher_exact_test(df_counts, time_threshold=18, testtype='less'):
    fisher_output = {}
    for primer in list(set(df_counts['primer'])):
        fisher_output[primer] = fisher_exact_test_by_primer(df_counts[df_counts['primer'] == primer],
                                                            time_threshold=time_threshold,
                                                            testtype=testtype)
    return fisher_output

def fisher_exact_test_replicate(df_counts, testtype='two-sided'):
    fisher_output = {}
    for primer in df_counts:
        fisher_output[primer] = fisher_exact_test_by_primer_replicate(df_counts[primer],
                                                                      testtype=testtype)
    return fisher_output

def fisher_exact_test_by_primer_replicate(df_counts_primer,testtype='two-sided'):
    df = df_counts_primer.sort_index()
    col_list= list(df)
    times = list(set([c[0] for c in col_list]))
    totals = [sum(df[c]) for c in col_list]
    for i,c in enumerate(col_list):
        df[(c[0],c[1]+0.1)] = totals[i] - df[c]

    #  Dict of how many replicates at each timepoint
    rep_dict = {}
    for key in list(df):
        if key[0] not in rep_dict:
            rep_dict[key[0]] = []
        if key[1] % 1 == 0:
            rep_dict[key[0]].append(key[1])
    for ti in single_rep_ti:
        del rep_dict[ti]
    for idx,row in df.iterrows():
        for ti in rep_dict:
            for index, rep_num1 in enumerate(rep_dict[ti]):
                for rep_num2 in rep_dict[ti][index+1:]:
                    oddsratio,pvalue = stats.fisher_exact([[row[(ti,rep_num1)],row[(ti,rep_num1)]],
                                    [row[(ti,rep_num2)],row[(ti,rep_num2)]]],
                                   alternative='two-sided')
                    key = '('+str(rep_num1)+','+str(rep_num2)+')'
                    df.at[idx,(ti,'oddsratio'+key)] = oddsratio
                    df.at[idx,(ti,'pvalue'+key)] = pvalue
    return df

def fisher_exact_test_by_primer(df_counts_primer,time_threshold=18,testtype='less'):
    df = df_counts_primer.sort_index()
    col_list= list(df)
    df['total'] = df[col_list].sum(axis=1)
    early_ti = [ti for ti in col_list if ti < time_threshold]
    late_ti = [ti for ti in col_list if ti >= time_threshold]
    df['early'] = df[early_ti].sum(axis=1)
    df['late'] = df[late_ti].sum(axis=1)
    df['late/early'] = df['late']/df['early']
    total_early = sum(df['early'])
    total_late = sum(df['late'])
    df['other_early'] = total_early - df['early']
    df['other_late'] = total_late - df['late']


    for idx,row in df.iterrows():
        oddsratio,pvalue = stats.fisher_exact([[row['early'],row['other_early']],
                    [row['late'],row['other_late']]],
                   alternative='less')
        df.at[idx,'oddsratio'] = oddsratio
        df.at[idx,'pvalue'] = pvalue
    df['oddsratio'] = np.reciprocal(df['oddsratio'])
    return df

def main():
    import argparse

    parser = argparse.ArgumentParser(
        description='Obtain counts from lineages to use in R expansion analysis.')
    parser.add_argument('--lineages', type=str,
                        help='path to json lineages file')
    parser.add_argument('--productive', action='store_false', default=True,
                        help='get counts only for productive lineages. (default: True)')
    parser.add_argument('--mergedcdr3', action='store_false', default=True,
                        help='use common cdr3 after merging (default: False)')
    parser.add_argument('--outfile', type=str,
                        help='path and name to save .csv file')
    args = parser.parse_args()

    if args.lineages is None:
        print('No input file')
        return
    if args.productive == True:
        print('Processing only productive lineages')
    else:
        print('Processing productive and unproductive lineages.')

    lineages = json_open(args.lineages)
    #  Get counts for productive lineages.
    csv_dict = create_csv_for_analysis(lineages['productive'], productive=args.productive,
                                       mergedcdr3=args.mergedcdr3)
    pd.DataFrame(csv_dict).to_csv(args.outfile, index=False)
    #make_csv_for_analysis(args.outfile, csv_dict)

    #  Get counts for unproductive lineages.
    csv_dict = create_csv_for_analysis(lineages['unproductive'], productive=False,
                                       mergedcdr3=args.mergedcdr3)
    pd.DataFrame(csv_dict).to_csv(args.outfile.replace('.csv', '_unproductive.csv'), index=False)
    #make_csv_for_analysis(args.outfile.replace('.csv', '_unproductive.csv'), csv_dict)

if __name__ == '__main__':
    main()
