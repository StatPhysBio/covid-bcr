import csv
from utils import *

def get_expansion_counts(lineages, patient_key, in_time_threshold, include_singletons=False):
    time_threshold = in_time_threshold
    primer_split_expansion = {"unique": {}, "abundance": {}}

    for key in lineages:
        lineage = lineages[key]
        uids = [ann['unique_ids'][0] for ann in lineage]

        lineage_primer = Counter([get_primer(u) for u in uids]).most_common(1)[0][0]
        if lineage_primer not in primer_split_expansion["unique"]:
            primer_split_expansion["unique"][lineage_primer] = {}
            primer_split_expansion["abundance"][lineage_primer] = {}


        unique_count = len(uids)
        abundance_counts = sum([get_dupcounts(u) for u in uids])
        times = np.array(CONST_DATA_DICT[int(patient_key)]['times']).astype(int)

        count_threshold = 1
        if include_singletons:
            count_threshold = 0

        for count_type in primer_split_expansion:
            counts = {}
            for ti in times:
                counts[ti] = sum([get_dupcounts(u) for u in uids
                              if get_time(u) == ti
                              and get_dupcounts(u) != count_threshold])
            primer_split_expansion[count_type][lineage_primer][key] = counts

    return primer_split_expansion

def create_data_csv(out_dir, data, in_patient, key1, key2):
    save_name = 'patient-' + in_patient + "_counts-" + key1 + "_primer-" + key2 + ".csv"
    csv_file = join(out_dir, save_name)
    col_titles = 'v_gene,j_gene,cdr3_length,'
    times = np.array(CONST_DATA_DICT[int(in_patient)]['times']).astype(int)
    for ti in times:
        col_titles+='clone_abundance_time_'+str(ti)+","
    col_titles = col_titles[:-1]+"\n"
    if len(times) == 3:
        with open(csv_file, 'w') as f:
            f.write(col_titles)
            for i,k in enumerate(data):
                time_data = []
                for ti in times:
                    time_data.append(data[k][ti])
                f.write("%s,%s,%s,%d,%d,%d\n"%(k[0],k[1],k[2],time_data[0],time_data[1],time_data[2]))
        print("3Wrote",csv_file)
    elif len(times) == 2:
        with open(csv_file, 'w') as f:
            f.write("v_gene,j_gene,cdr3_length,clone_early_abundance,clone_late_abundance\n")
            for i,k in enumerate(data):
                time_data = []
                for ti in times:
                    time_data.append(data[k][ti])
                f.write("%s,%s,%s,%d,%d\n"%(k[0],k[1],k[2],time_data[0],time_data[1]))
        print("2Wrote",csv_file)

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
