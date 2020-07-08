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
                             get_lineage_progenitor_cdr3)
from utils import *

def get_counts_by_time(headers: list, headers_unique_counts: list, times: list) -> dict:
    """Obtains unique, abundance, and singleton counts per time in a lineage.

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
        at all times in a lineage.
    lineage_primer : str
        The primer used by the majority of sequences in a lineage. (They should
        all have the same primer if the pipeline, which was supposed to remove
        sequences with mismatched primers, was run correctly.)
    """

    count_dict = {'unique': {},
                  'unique_merged': {},
                  'abundance': {},
                  'singleton': {}}
    vprimers = [None]*len(headers_unique_counts)

    for ti in times:
        count_dict['unique'][ti] = 0
        count_dict['unique_merged'][ti] = 0
        count_dict['abundance'][ti] = 0
        count_dict['singleton'][ti] = 0

    for i,u in enumerate(headers_unique_counts):
        vprimers[i] = get_vprimer(u)
        ti = get_time(u)
        count_dict['unique_merged'][ti] += 1

    lineage_primer = Counter(vprimers).most_common(1)[0][0]

    for u in headers:
        ti = get_time(u)
        count_dict['unique'][ti] += 1
        abundance = get_abundance(u)
        count_dict['abundance'][ti] += abundance
        if abundance == 1:
            count_dict['singleton'][ti] += 1

    return count_dict,lineage_primer

def get_counts_by_time_replicate(headers: list, headers_unique_counts: list, times: list) -> dict:
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

    count_dict = {'unique': {},
                  'unique_merged': {},
                  'abundance': {},
                  'singleton': {}}
    vprimers = [None]*len(headers_unique_counts)

    for ti in times:
        #  Instead of running a loop over all the data and obtaining
        #  the replicates, we will put 0's for all possible replicates. The
        #  replicates which don't exist will be removed quickly in another function.
        for r in range(0,5):
            count_dict['unique'][(ti,r)] = 0
            count_dict['unique_merged'][(ti,r)] = 0
            count_dict['abundance'][(ti,r)] = 0
            count_dict['singleton'][(ti,r)] = 0

    for i,u in enumerate(headers_unique_counts):
        vprimers[i] = get_vprimer(u)
        ti = get_time(u)
        rep = get_replicate(u)
        count_dict['unique_merged'][(ti,rep)] += 1

    lineage_primer = Counter(vprimers).most_common(1)[0][0]

    for u in headers:
        ti = get_time(u)
        rep = get_replicate(u)
        count_dict['unique'][(ti,rep)] += 1
        abundance = get_abundance(u)
        count_dict['abundance'][(ti,rep)] += abundance
        if abundance == 1:
            count_dict['singleton'][(ti,rep)] += 1

    return count_dict,lineage_primer

def get_counts(in_lineages: dict, productive: bool = True,
               replicate: bool = False, abstar: bool = True) -> pd.DataFrame:
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

    primer_split_counts = {'unique': {},
                           'unique_merged': {},
                           'abundance': {},
                           'singleton': {}}
    if replicate:
        count_func = get_counts_by_time_replicate
    else:
        count_func = get_counts_by_time

    lineages = denest_lineages(in_lineages, productive=productive)
    patient_key = get_patient(lineages[list(lineages.keys())[0]][0]['seq_id'])

    #  Looks up the times from the CONST_DATA_DICT in utils.
    try:
        times = CONST_DATA_DICT[str(patient_key)]['sample day']
    except:
        print('Patient key given is not present.')
        return

    for key in lineages:
        lineage = lineages[key]
        lineage_unique_counts = merge_replicates_within_lineage(lineage)

        #  Annotations come from abstar.
        if abstar:
            headers = [ann['seq_id'] for ann in lineage]
            headers_unique_counts = [ann['seq_id'] for ann in lineage_unique_counts]
        #  Annotations come from partis.
        else:
            headers = [ann['unique_ids'][0] for ann in lineage]
            headers_unique_counts = [ann['unique_ids'][0] for ann in lineage_unique_counts]

        count_dict,lineage_primer = count_func(headers, headers_unique_counts, times)
        if lineage_primer not in primer_split_counts['unique']:
            primer_split_counts['unique'][lineage_primer] = {}
            primer_split_counts['unique_merged'][lineage_primer] = {}
            primer_split_counts['abundance'][lineage_primer] = {}
            primer_split_counts['singleton'][lineage_primer] = {}
        for count_type in primer_split_counts:
            primer_split_counts[count_type][lineage_primer][key] = count_dict[count_type]

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
        df_counts[count_type] = df_counts[count_type].loc[:, (df_counts[count_type] != 0).any(axis=0)]
        #  Resetting indices causes the primer, v gene, j gene, cdr3 length, and cluster id
        #  to become columns. This is much more useful when using the fisher exact test
        #  for isolating the counts by primer.
        df_counts[count_type] = df_counts[count_type].reset_index()
        #  Rename the columns so that they are meaningful.
        df_counts[count_type] = df_counts[count_type].rename(columns={'level_0': 'primer',
                                                                      'level_1': 'v_gene',
                                                                      'level_2': 'j_gene',
                                                                      'level_3': 'cdr3_length',
                                                                      'level_4': 'cluster_id'})
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

def create_csv_for_analysis(lineages: dict, productive: bool = True, replicate: bool = True) -> dict:
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

    csv_dict = {'patient':[], 'v_gene': [], 'j_gene': [], 'cdr3_length': [],
                'common_cdr3':[], 'cluster_id':[], 'time':[], 'replicate':[],
                'abundance':[], 'unique': [], 'unique_merged': [], 'singleton': [],
                'productive':[],'primer':[]}

    df_counts = get_counts(lineages, productive=productive,
                           replicate=replicate, abstar=True)
    ti_for_lins = get_columns_with_ints(df_counts['abundance'])
    lins = df_counts['abundance'][['v_gene','j_gene','cdr3_length','cluster_id','primer']].values.tolist()
    patient = get_patient(lineages[lins[0][0]][lins[0][1]][lins[0][2]][lins[0][3]][0]['seq_id'])

    common_cdr3 = {}
    for lin in lins:
        lineage = lineages[lin[0]][lin[1]][lin[2]][lin[3]]
        common_cdr3[tuple(lin)] = get_lineage_progenitor_cdr3(lineage, nt=False)
    print('# lins', len(lins))
    print('Lins with unproductive CDR3s',sum([1 for key in common_cdr3 if common_cdr3[key] == '']))

    for i, lin in enumerate(lins):
        #  Don't save any productive lineages with unproductive progenitors.
        if common_cdr3[tuple(lin)] == '' and productive==True:
            continue
        for ti in ti_for_lins:
            csv_dict['primer'].append(lin[4])
            csv_dict['patient'].append(patient)
            csv_dict['v_gene'].append(lin[0])
            csv_dict['j_gene'].append(lin[1])
            csv_dict['cdr3_length'].append(int(lin[2]))
            csv_dict['common_cdr3'].append(common_cdr3[tuple(lin)])
            csv_dict['cluster_id'].append(int(lin[3]))
            csv_dict['time'].append(ti[0])
            csv_dict['replicate'].append(ti[1])
            csv_dict['abundance'].append(df_counts['abundance'].iloc[i][ti])
            csv_dict['unique'].append(df_counts['unique'].iloc[i][ti])
            csv_dict['unique_merged'].append(df_counts['unique_merged'].iloc[i][ti])
            csv_dict['singleton'].append(df_counts['singleton'].iloc[i][ti])
            csv_dict['productive'].append(productive)
    return csv_dict

def make_csv_for_analysis(savename: str, csv_dict: dict, replicate: bool = True) -> None:
    """Saves the necessary information to a csv for R expansion analysis.

    Parameters
    ----------
    savename : str
        Path and name of csv file to be saved.
    csv_dict : dict
        Dictionary of data necessary for expansion analysis.
    replicate : bool, optional
        Bool to specify whether or not the replicate column is necessary.

    Returns
    -------
    None
    """

    col_list = ['patient','primer','productive',
                'v_gene','j_gene','cdr3_length','cluster_id',
                'common_cdr3','time','replicate',
                'abundance','unique','unique_merged','singleton']

    if not replicate:
        col_list.remove('replicate')

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
    #  Delete any timepoints with only one rep
    single_rep_ti = [ti for ti in rep_dict if len(rep_dict[ti]) < 2]
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
    parser.add_argument('--productive', action='store_true', default=True,
                        help='get counts only for productive lineages. (default: True)')
    parser.add_argument('--replicate', action='store_true', default=True,
                        help='get counts separated by replicate (default: True)')
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
    if args.replicate == True:
        print('Processing replicates separately')
    else:
        print('Processing replicates jointly.')

    lineages = json_open(args.lineages)
    #  Get counts for productive lineages.
    csv_dict = create_csv_for_analysis(lineages['productive'], productive=args.productive,
                                       replicate=args.replicate)
    make_csv_for_analysis(args.outfile, csv_dict, replicate=args.replicate)

    #  Get counts for unproductive lineages.
    csv_dict = create_csv_for_analysis(lineages['unproductive'], productive=False,
                                       replicate=args.replicate)
    make_csv_for_analysis(args.outfile.replace('.csv', '_unproductive.csv'),
                          csv_dict, replicate=args.replicate)

if __name__ == '__main__':
    main()
