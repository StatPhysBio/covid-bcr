import csv
from utils import *
import scipy.stats as stats

def get_expansion_counts(lineages, patient_key, include_singletons=False, slc=True):
    primer_split_expansion = {"unique": {},
                              "abundance": {},
                              "singleton": {}}

    if not slc:
        for key in lineages:
            lineage = lineages[key]
            uids = [ann['unique_ids'][0] for ann in lineage]

            lineage_primer = Counter([get_vprimer(u) for u in uids]).most_common(1)[0][0]
            if lineage_primer not in primer_split_expansion["unique"]:
                primer_split_expansion["unique"][lineage_primer] = {}
                primer_split_expansion["abundance"][lineage_primer] = {}
                primer_split_expansion["singleton"][lineage_primer] = {}

            unique_count = len(uids)
            abundance_counts = sum([get_abundance(u) for u in uids])
            times = np.array(CONST_DATA_DICT[int(patient_key)]['times']).astype(int)

            count_threshold = 1
            if include_singletons:
                count_threshold = 0

            for count_type in primer_split_expansion:
                if count_type == 'unique':
                    counts = {}
                    for ti in times:
                        counts[ti] = sum([1 for u in uids
                                      if get_time(u) == ti
                                      and get_abundance(u) != count_threshold])
                    primer_split_expansion[count_type][lineage_primer][key] = counts
                elif count_type == 'abundance':
                    counts = {}
                    for ti in times:
                        counts[ti] = sum([get_abundance(u) for u in uids
                                      if get_time(u) == ti
                                      and get_abundance(u) != count_threshold])
                    primer_split_expansion[count_type][lineage_primer][key] = counts
                elif count_type == 'singleton':
                    counts = {}
                    for ti in times:
                        counts[ti] = sum([get_abundance(u) for u in uids
                                      if get_time(u) == ti
                                      and get_abundance(u) == 1])
                    primer_split_expansion[count_type][lineage_primer][key] = counts

    else:
        for key in lineages:
            for rank in lineages[key]:
                lineage = lineages[key][rank]
                uids = [ann['unique_ids'][0] for ann in lineage]

                lineage_primer = Counter([get_vprimer(u) for u in uids]).most_common(1)[0][0]
                if lineage_primer not in primer_split_expansion["unique"]:
                    primer_split_expansion["unique"][lineage_primer] = {}
                    primer_split_expansion["abundance"][lineage_primer] = {}
                    primer_split_expansion["singleton"][lineage_primer] = {}

                if key not in primer_split_expansion["unique"][lineage_primer]:
                    primer_split_expansion["unique"][lineage_primer][key] = {}
                    primer_split_expansion["abundance"][lineage_primer][key] = {}
                    primer_split_expansion["singleton"][lineage_primer][key] = {}

                unique_count = len(uids)
                abundance_counts = sum([get_abundance(u) for u in uids])
                times = np.array(CONST_DATA_DICT[int(patient_key)]['times']).astype(int)

                count_threshold = 1
                if include_singletons:
                    count_threshold = 0

                for count_type in primer_split_expansion:
                    if count_type == 'unique':
                        counts = {}
                        for ti in times:
                            counts[ti] = sum([1 for u in uids
                                          if get_time(u) == ti
                                          and get_abundance(u) != count_threshold])
                        primer_split_expansion[count_type][lineage_primer][key][rank] = counts
                    elif count_type == 'singleton':
                        counts = {}
                        for ti in times:
                            counts[ti] = sum([get_abundance(u) for u in uids
                                          if get_time(u) == ti
                                          and get_abundance(u) == 1])
                        primer_split_expansion[count_type][lineage_primer][key][rank] = counts
                    elif count_type == 'abundance':
                        counts = {}
                        for ti in times:
                            counts[ti] = sum([get_abundance(u) for u in uids
                                          if get_time(u) == ti
                                          and get_abundance(u) != count_threshold])
                        primer_split_expansion[count_type][lineage_primer][key][rank] = counts

    return primer_split_expansion

def fisher_exact_test(in_expansion, slc=False, time_threshold=18, testtype='less'):
    fisher_output = {}
    for primer in in_expansion:
        fisher_output[primer] = fisher_exact_test_by_primer(in_expansion[primer], slc=slc,
                                                         time_threshold=time_threshold,testtype=testtype)
    return fisher_output

def fisher_exact_test_by_primer(primer_bin, slc=False, time_threshold=18,testtype='less'):
    #  Get counts
    counts = {}
    if not slc:
        for key in primer_bin:
            for ti in primer_bin[key]:
                if int(ti) not in counts:
                    counts[int(ti)] = []
                counts[int(ti)].append(primer_bin[key][ti])
    else:
        for key in primer_bin:
            for rank in primer_bin[key]:
                for ti in primer_bin[key][rank]:
                    if int(ti) not in counts:
                        counts[int(ti)] = []
                    counts[int(ti)].append(primer_bin[key][rank][ti])

    #  Split into late and early
    fisher_input = {"early":[],"late":[], "other early": [], "other late":[]}
    for i,c in enumerate(counts[list(counts.keys())[0]]):
        fisher_input["early"].append(sum([counts[ti][i] for ti in counts
                                          if ti < time_threshold]))
        fisher_input["late"].append(sum([counts[ti][i] for ti in counts
                                          if ti >= time_threshold]))

    #  Perform fisher exact test
    fisher_results = {"oddsratio":[],"pvalue":[]}
    for i,e_value in enumerate(fisher_input['early']):
        l_value = fisher_input['late'][i]
        other_early = sum([ecount
                           for k,ecount in enumerate(fisher_input['early'])
                           if k != i])
        fisher_input["other early"].append(other_early)
        other_late = sum([lcount
                           for k,lcount in enumerate(fisher_input['late'])
                           if k != i])
        fisher_input["other late"].append(other_late)
        oddsratio, pvalue = stats.fisher_exact([[e_value, other_early],
                                                [l_value, other_late]],
                                               alternative=testtype)
        fisher_results['pvalue'].append(pvalue)
        fisher_results['oddsratio'].append(oddsratio)

    #  Set 0 to be extremely small number
    #  so it shows up on log pvalue plot
    #  fisher_results['pvalue'] = np.array(fisher_results['pvalue'])
    #  fisher_results['pvalue'][fisher_results['pvalue'] == 0.0] = 1e-340
    fisher_results['early'] = np.array(fisher_input['early'])
    fisher_results['late'] = np.array(fisher_input['late'])
    fisher_results['other early'] = np.array(fisher_input['other early'])
    fisher_results['other late'] = np.array(fisher_input['other late'])
    fisher_results['pvalue'] = np.array(fisher_results['pvalue'])
    fisher_results['oddsratio'] = np.reciprocal(np.array(fisher_results['oddsratio']))

    if not slc:
        fisher_results['vjl'] = list(primer_bin.keys())
    else:
        fisher_results['vjl'] = [vjl
                                for vjl in primer_bin.keys()
                                for rank in primer_bin[vjl]]
        fisher_results['vjl+rank'] = [list(vjl)+[rank]
                                      for vjl in primer_bin.keys()
                                      for rank in primer_bin[vjl]]
    return fisher_results

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
