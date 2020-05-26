import csv
from utils import *
import scipy.stats as stats
import pandas as pd

def get_counts_by_time(uids, times):
    count_dict = {"unique": {},
                 "abundance": {},
                 "singleton": {}}
    vprimers = [None]*len(uids)
    for ti in times:
        count_dict['unique'][ti] = 0
        count_dict['abundance'][ti] = 0
        count_dict['singleton'][ti] = 0

    for i,u in enumerate(uids):
        vprimers[i] = get_vprimer(u)
        ti = get_time(u)
        count_dict['unique'][ti] += 1
        abundance = get_abundance(u)
        count_dict['abundance'][ti] += abundance
        if abundance == 1:
            count_dict['singleton'][ti] += 1
    lineage_primer = Counter(vprimers).most_common(1)[0][0]
    return count_dict,lineage_primer

def get_counts_by_time_replicate(uids, times):
    count_dict = {"unique": {},
                 "abundance": {},
                 "singleton": {}}
    vprimers = [None]*len(uids)
    for ti in times:
        for r in range(0,3):
            count_dict['unique'][(ti,r)] = 0
            count_dict['abundance'][(ti,r)] = 0
            count_dict['singleton'][(ti,r)] = 0

    for i,u in enumerate(uids):
        vprimers[i] = get_vprimer(u)
        ti = get_time(u)
        rep = get_replicate(u)
        count_dict['unique'][(ti,rep)] += 1
        abundance = get_abundance(u)
        count_dict['abundance'][(ti,rep)] += abundance
        if abundance == 1:
            count_dict['singleton'][(ti,rep)] += 1
    lineage_primer = Counter(vprimers).most_common(1)[0][0]
    return count_dict,lineage_primer

def get_primer_split_counts(lineages, patient_key, rep=False):
    primer_split_counts = {"unique": {},
                           "abundance": {},
                           "singleton": {}}
    if rep:
        count_func = get_counts_by_time_replicate
    else:
        count_func = get_counts_by_time
    if type(patient_key) == list:
        times = []
        for pkey in patient_key:
            times += CONST_DATA_DICT[str(pkey)]['sample day']
    else:
        times = CONST_DATA_DICT[str(patient_key)]['sample day']

    for key in lineages:
        lineage = lineages[key]
        uids = [ann['unique_ids'][0] for ann in lineage]
        count_dict,lineage_primer = count_func(uids,times)
        if lineage_primer not in primer_split_counts["unique"]:
            primer_split_counts["unique"][lineage_primer] = {}
            primer_split_counts["abundance"][lineage_primer] = {}
            primer_split_counts["singleton"][lineage_primer] = {}
        for count_type in primer_split_counts:
            primer_split_counts[count_type][lineage_primer][key] = count_dict[count_type]

    #  Convert to dataframe
    df_counts = {}
    for count_type in primer_split_counts:
        df_counts[count_type] = {}
        for lineage_primer in primer_split_counts[count_type]:
            df_counts[count_type][lineage_primer] = pd.DataFrame.from_dict(primer_split_counts[count_type][lineage_primer],
                                                                                      orient='index',dtype=np.uint32)
            df = df_counts[count_type][lineage_primer]
            df = df.loc[:, (df != 0).any(axis=0)]
            df_counts[count_type][lineage_primer] = df
    return df_counts


def fisher_exact_test(df_counts, time_threshold=18, testtype='less'):
    fisher_output = {}
    for primer in df_counts:
        fisher_output[primer] = fisher_exact_test_by_primer(df_counts[primer],
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
                    key = "("+str(rep_num1)+","+str(rep_num2)+")"
                    df.at[idx,(ti,'oddsratio'+key)] = oddsratio
                    df.at[idx,(ti,'pvalue'+key)] = pvalue
    #for ti in times:
    #    col_ti = list(df[ti])
    #    for c in col_ti:
    #        if 'oddsratio' in c:
    #            df[(ti,c)] = np.reciprocal(df[(ti,c)])
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

def write_lineage_info(fet,patient,rank,linkey,savedir):
    col_titles = 'v_gene,j_gene,cdr3_length,pvalue,late_counts/early_counts,oddsratio\n'
    v=linkey[0]
    j=linkey[1]
    cl=linkey[2]
    pvalue=fet.at[linkey,'pvalue']
    fold_expansion=fet.at[linkey,'late/early']
    oddsratio=fet.at[linkey,'oddsratio']
    with open(savedir+'patient-'+patient+'_lineage-'+str(rank)+'.csv','w') as f:
        f.write(col_titles)
        f.write("%s,%s,%s,%.16e,%.16e,%.16e,"%(v,j,cl,
                                               pvalue,
                                               fold_expansion,
                                               oddsratio))
def get_expanded_lineages(df_fet, pvalue_thresh=1e-200):
    expanded_lineages = []
    for primer in df_fet:
        conditions = ((df_fet[primer]['pvalue'] < pvalue_thresh)
                      & (df_fet[primer]['early'] != 0)
                      & (df_fet[primer]['late'] != 0))
        expanded_linkeys = df_fet[primer][conditions].index.values.tolist()
        print(expanded_linkeys)
        expanded_linkeys = [linkey for linkey in expanded_linkeys]
        expanded_lineages += expanded_linkeys
    return expanded_lineages

def make_trees_for_expanded_lineages(lineages,patient,savedir):
    from get_input_for_trees import make_fastas_for_trees
    df_counts = get_primer_split_counts(lineages, patient)['abundance']
    df_fet = fisher_exact_test(df_counts)
    expanded_lineages = []
    sizes = []
    for primer in df_fet:
        conditions = ((df_fet[primer]['pvalue'] < 1e-200)
                      & (df_fet[primer]['early'] != 0)
                      & (df_fet[primer]['late'] != 0))
        expanded_linkeys = df_fet[primer][conditions].index.values.tolist()
        expanded_linkeys = [tuple(list(linkey) + [primer]) for linkey in expanded_linkeys]
        expanded_lineages += expanded_linkeys
        sizes += df_fet[primer][conditions]['total'].tolist()

    sizes.sort(reverse=True)
    for linkey in expanded_lineages:
        lin = lineages[(linkey[0:3])][linkey[3]]
        sizerank = sizes.index(df_fet[linkey[4]].at[linkey[0:4],'total'])
        make_fastas_for_trees(lin, sizerank, patient,
                              True, False,savedir)
        write_lineage_info(df_fet[linkey[4]],patient,sizerank,linkey[0:4],savedir)

#  For comparison to Oxford database
def make_csv_file(savename, in_fet_output,in_expansion,slc=False):
    print(in_expansion.keys())
    if not slc:
        tkey = list(in_expansion['singleton']['IGHV1-F'].keys())[0]
        times = list(in_expansion['singleton']['IGHV1-F'][tkey].keys())
    else:
        tvjl = list(in_expansion['singleton']['IGHV1-F'].keys())[0]
        trank = list(in_expansion['singleton']['IGHV1-F'][tvjl].keys())[0]
        times = list(in_expansion['singleton']['IGHV1-F'][tvjl][trank].keys())
    ti_col_names = ['_unique',"_abundance","_singletons"]
    ti_cols = []
    if len(times) != 3:
        times.append("")
    for s in ti_col_names:
        for ti in times:
            ti_cols.append(str(ti)+s)
    col_end = ",".join(ti_cols)
    col_titles = 'v_gene,j_gene,cdr3_length,pvalue,late_counts/early_counts,oddsratio,'+col_end+"\n"

    fet_output = sort_dict(in_fet_output)
    with open(savename,'w') as f:
        f.write(col_titles)
        for primer in fet_output:
            for index,lin in enumerate(fet_output[primer]['vjl']):
                fold_expansion = fet_output[primer]['fet']['late'][index]/fet_output[primer]['fet']['early'][index]
                oddsratio = fet_output[primer]['fet']['oddsratio'][index]
                j, cl, linrank = '','',''
                if len(lin) > 4:
                    v = lin
                elif len(lin) == 2:
                    v = lin[0]
                    j = lin[1]
                elif len(lin) == 3:
                    v = lin[0]
                    j = lin[1]
                    cl = lin[2]
                elif len(lin) == 4:
                    v = lin[0]
                    j = lin[1]
                    cl = lin[2]
                f.write("%s,%s,%s,%.16e,%.16e,%.16e,"%(v,j,cl,
                                                        fet_output[primer]['fet']['pvalue'][index],
                                                        fold_expansion,
                                                        oddsratio))
                if not slc:
                    c = 0
                    for count_type in in_expansion:
                        for ti in times:
                            if not ti=="":
                                f.write("%i"%in_expansion[count_type][primer][lin][ti])
                            c+=1
                            if c != len(in_expansion)*len(times):
                                f.write(",")
                    f.write("\n")
                else:
                    linrank = fet_output[primer]['vjl+rank'][index][3]
                    c = 0
                    for count_type in in_expansion:
                        for ti in times:
                            if not ti=="":
                                f.write("%i"%in_expansion[count_type][primer][lin][linrank][ti])
                            c+=1
                            if c != len(in_expansion)*len(times):
                                f.write(",")
                    f.write("\n")
    print("Finished making", savename)
#  Deprecated.
#def create_data_csv(out_dir, data, in_patient, key1, key2, extra_info=""):
#    save_name = 'patient-' + in_patient + "_counts-" + key1 + "_primer-" + key2 + "_"+extra_info+".csv"
#    csv_file = join(out_dir, save_name)
#    col_titles = 'v_gene,j_gene,cdr3_length,'
#    times = np.array(CONST_DATA_DICT[int(in_patient)]['times']).astype(int)
#    for ti in times:
#        col_titles+='clone_abundance_time_'+str(ti)+","
#    col_titles = col_titles[:-1]+"\n"
#    if len(times) == 3:
#        with open(csv_file, 'w') as f:
#            f.write(col_titles)
#            for i,k in enumerate(data):
#                time_data = []
#                for ti in times:
#                    time_data.append(data[k][ti])
#                f.write("%s,%s,%s,%d,%d,%d\n"%(k[0],k[1],k[2],time_data[0],time_data[1],time_data[2]))
#        print("3Wrote",csv_file)
#    elif len(times) == 2:
#        with open(csv_file, 'w') as f:
#            f.write("v_gene,j_gene,cdr3_length,clone_early_abundance,clone_late_abundance\n")
#            for i,k in enumerate(data):
#                time_data = []
#                for ti in times:
#                    time_data.append(data[k][ti])
#                f.write("%s,%s,%s,%d,%d\n"%(k[0],k[1],k[2],time_data[0],time_data[1]))
#        print("2Wrote",csv_file)

#def main():
#    import argparse
#    parser = argparse.ArgumentParser(
#        description='Create csv file for fisher exact test for expansion')
#    parser.add_argument('--lineages', type=str, help='path to lineage file')
#    parser.add_argument('--outdir', type=str, help='path to outfile')
#    parser.add_argument('--patient', type=str, help='patient num')
#    parser.add_argument('--time_threshold', type=int, help='time used to split late and early')
#    parser.add_argument('--singletons', action='store_true', help='include singletons counts')
#    parser.add_argument('--unique', action='store_true', help='use unique counts')
#    parser.add_argument('--abundance', action='store_true', help='use abundance counts')
#    args = parser.parse_args()
#
#    lineages = unpickle(args.lineages)
#    expansion_counts = get_expansion_counts(lineages, args.patient, args.time_threshold, args.singletons)
#
#    if args.unique:
#        for primer in expansion_counts['unique']:
#            create_data_csv(args.outdir, expansion_counts['abundance'][primer],args.patient,'abundance',primer)
#
#    if args.abundance:
#        for primer in expansion_counts['abundance']:
#            create_data_csv(args.outdir, expansion_counts['abundance'][primer],args.patient,'abundance',primer)
#
#if __name__ == '__main__':
#    main()
