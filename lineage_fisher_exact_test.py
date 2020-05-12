import csv
from utils import *
import scipy.stats as stats

def get_expansion_counts(lineages, patient_key, include_singletons=False, slc=True):
    primer_split_expansion = {"unique": {}, "abundance": {}}

    if not slc:
        for key in lineages:
            lineage = lineages[key]
            uids = [ann['unique_ids'][0] for ann in lineage]

            lineage_primer = Counter([get_vprimer(u) for u in uids]).most_common(1)[0][0]
            if lineage_primer not in primer_split_expansion["unique"]:
                primer_split_expansion["unique"][lineage_primer] = {}
                primer_split_expansion["abundance"][lineage_primer] = {}


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
                else:
                    counts = {}
                    for ti in times:
                        counts[ti] = sum([get_abundance(u) for u in uids
                                      if get_time(u) == ti
                                      and get_abundance(u) != count_threshold])
                    primer_split_expansion[count_type][lineage_primer][key] = counts
    else:
        for key in lineages:
            print(key)
            for rank in lineages[key]:
                print(rank)
                lineage = lineages[key][rank]
                uids = [ann['unique_ids'][0] for ann in lineage]

                lineage_primer = Counter([get_vprimer(u) for u in uids]).most_common(1)[0][0]
                if lineage_primer not in primer_split_expansion["unique"]:
                    primer_split_expansion["unique"][lineage_primer] = {}
                    primer_split_expansion["abundance"][lineage_primer] = {}

                if key not in primer_split_expansion["unique"][lineage_primer]:
                    primer_split_expansion["unique"][lineage_primer][key] = {}
                    primer_split_expansion["abundance"][lineage_primer][key] = {}

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
                    else:
                        counts = {}
                        for ti in times:
                            counts[ti] = sum([get_abundance(u) for u in uids
                                          if get_time(u) == ti
                                          and get_abundance(u) != count_threshold])
                        primer_split_expansion[count_type][lineage_primer][key][rank] = counts

    return primer_split_expansion

def fisher_exact_test(primer_bin, slc=False, time_threshold=18,testtype='two-sided'):
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
                print(key, rank)
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
#    fisher_results['pvalue'] = np.array(fisher_results['pvalue'])
#    fisher_results['pvalue'][fisher_results['pvalue'] == 0.0] = 1e-300

    return (np.array(fisher_input['early']),
            np.array(fisher_input['late']),
            np.array(fisher_input['other early']),
            np.array(fisher_input['other late']),
            np.array(fisher_results['pvalue']),
            np.reciprocal(np.array(fisher_results['oddsratio'])))

def get_vjls_and_fet(lineages, patient, timecutoff,slc=False):
    print(slc)
    expansion = get_expansion_counts(lineages, patient, include_singletons=True, slc=slc)['abundance']
    fisher_output = {}

    for primer in expansion:
        fisher_output[primer] = {}
        fisher_output[primer]['fet'] =  fisher_exact_test(expansion[primer],
                                                          time_threshold=timecutoff,
                                                          testtype='less',slc=slc)
        if not slc:
            fisher_output[primer]['vjl'] = list(expansion[primer].keys())
        else:
            fisher_output[primer]['vjl'] = [key
                                            for key in expansion[primer].keys()
                                            for item in expansion[primer][key]]
            fisher_output[primer]['vjl+rank'] = [key
                                            for key in expansion[primer].keys()
                                            for item in expansion[primer][key]]
    return fisher_output

#  For comparison to Oxford database
def make_csv_file(savename, in_fet_output):
    col_titles = 'v_gene,j_gene,cdr3_length,pvalue,late_counts/early_counts,oddsratio\n'
    fet_output = sort_dict(in_fet_output)
    with open(savename,'w') as f:
        f.write(col_titles)
        for primer in fet_output:
            for index,lin in enumerate(fet_output[primer]['vjl']):
                if fet_output[primer]['fet'][0][index] == 0 or fet_output[primer]['fet'][1][index] == 0:
                    continue
                fold_expansion = fet_output[primer]['fet'][1][index]/fet_output[primer]['fet'][0][index]
                oddsratio = fet_output[primer]['fet'][5][index]
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
                f.write("%s,%s,%s,%.16e,%.16e,%.16e\n"%(v,j,cl,
                                                        fet_output[primer]['fet'][4][index],
                                                        fold_expansion,
                                                        oddsratio))

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

def main():
    import argparse
    parser = argparse.ArgumentParser(
        description='Create csv file for fisher exact test for expansion')
    parser.add_argument('--lineages', type=str, help='path to lineage file')
    parser.add_argument('--outdir', type=str, help='path to outfile')
    parser.add_argument('--patient', type=str, help='patient num')
    parser.add_argument('--time_threshold', type=int, help='time used to split late and early')
    parser.add_argument('--singletons', action='store_true', help='include singletons counts')
    parser.add_argument('--unique', action='store_true', help='use unique counts')
    parser.add_argument('--abundance', action='store_true', help='use abundance counts')
    args = parser.parse_args()

    lineages = unpickle(args.lineages)
    expansion_counts = get_expansion_counts(lineages, args.patient, args.time_threshold, args.singletons)

    if args.unique:
        for primer in expansion_counts['unique']:
            create_data_csv(args.outdir, expansion_counts['abundance'][primer],args.patient,'abundance',primer)

    if args.abundance:
        for primer in expansion_counts['abundance']:
            create_data_csv(args.outdir, expansion_counts['abundance'][primer],args.patient,'abundance',primer)

if __name__ == '__main__':
    main()
