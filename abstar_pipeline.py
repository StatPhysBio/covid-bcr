from utils import *
from error_correct import *
from vjl_slc import vjl_slc
from lineage_fisher_exact_test import *

def load_abstar(abstar_file):
    num_N = 0
    len_fail = 0
    d_fail = 0
    cdr3_fail = 0
    passed = 0
    annotations = {'productive':{},
                   'unproductive':{},
                   'stop_in_cdr3':{},
                   'other':{}}
    for line in open(abstar_file, 'r'):
        if len(line) < 10:
            len_fail += 1
            continue
        if 'd_gene' not in line:
            d_fail += 1
            continue
        data_dict = json.loads(line)
        if 'cdr3_nt' not in data_dict:
            cdr3_fail += 1
            continue
        if data_dict['raw_input'].strip("N").count("N") != 0:
            num_N += 1
            continue
        passed += 1
        seq_id = data_dict['seq_id'].split("|")[0]
        if data_dict['prod'] == 'yes':
            annotations['productive'][seq_id] = data_dict
        elif data_dict['junction_in_frame'] == 'no':
            annotations['unproductive'][seq_id] = data_dict
        elif "*" in data_dict['junc_aa']:
            annotations['stop_in_cdr3'][seq_id] = data_dict
        else:
            annotations['other'][seq_id] = data_dict
    print("Length_fail",len_fail,
          "\nD_gene_fail",d_fail,
          "\nN_filtered", num_N,
          "\nCDR3_fail",cdr3_fail,
          '\nabstar_pass', passed)
    return annotations

def shm_indels(in_ann):
    return (len(in_ann['vdj_nt']) != len(in_ann['gapped_vdj_nt'])
       or len(in_ann['vdj_germ_nt']) != len(in_ann['gapped_vdj_germ_nt']))

def abstar_cdr3(in_a):
    try:
        cdr3_start = in_a['vdj_nt'].index(in_a['cdr3_nt'])
        cdr3_end = cdr3_start + in_a['region_len_nt']['cdr3']
        return cdr3_start,cdr3_end
    except:
        return None

def run_error_correction(annotations, keys):
    if ('productive' not in keys and 'unproductive' not in keys) or ('productive' in keys and 'unproductive' in keys):
        return

    in_h = []
    in_seqs = []
    for key in keys:
        for seq_key in annotations[key]:
            in_h.append(annotations[key][seq_key]['seq_id'])
            in_seqs.append(annotations[key][seq_key]['raw_input'])

    before_ec_unique = len(in_h)
    before_ec_abundance = sum([get_abundance(u) for u in in_h])
    marg_uids = []

    groups = group_data(in_h,in_seqs)
    for key in groups:
        marg_out_uids, marg_out_seqs, marg_pc = error_correct_marginal(groups[key]['uids'],groups[key]['sequences'])
        marg_uids += marg_out_uids

        tot_out_uids, tot_out_seqs, total_pc = error_correct_total(marg_out_uids,marg_out_seqs,d_tol=2)
        groups[key]['uids'] = tot_out_uids
        groups[key]['sequences'] = tot_out_seqs

    final_uids = []
    final_seqs = []
    for key in groups:
        final_uids += groups[key]['uids']
        final_seqs += groups[key]['sequences']

    #  Report counts
    marg_unique = 0
    marg_abundance = 0
    for u in marg_uids:
        abun = get_abundance(u)
        if abun > 0:
            marg_unique += 1
            marg_abundance += abun
    final_unique = 0
    final_abundance = 0
    for u in final_uids:
        abun = get_abundance(u)
        if abun > 0:
            final_unique += 1
            final_abundance += abun

    updated_dict_unique = 0
    updated_dict_abundance = 0
    if 'productive' in keys:
        prefix = 'prod_'
        for uid in final_uids:
            seqkey = uid.split("|")[0]
            try:
                abun = get_abundance(uid)
                if abun > 0:
                    updated_dict_unique += 1
                    updated_dict_abundance += abun
                    annotations['productive'][seqkey]['seq_id'] = uid
                else:
                    del annotations['productive'][seqkey]
            except:
                pass

    else:
        prefix = 'unprod_'
        for uid in final_uids:
            seqkey = uid.split("|")[0]
            try:
                abun = get_abundance(uid)
                if abun > 0:
                    updated_dict_unique += 1
                    updated_dict_abundance += abun
                    annotations['unproductive'][seqkey]['seq_id'] = uid
                else:
                    del annotations['unproductive'][seqkey]
            except:
                pass

    print(prefix + 'before_ec', before_ec_unique, before_ec_abundance)
    print(prefix + 'after_marginal', marg_unique, marg_abundance)
    print(prefix + 'after_total', final_unique, final_abundance)
    print(prefix + 'after_type_check', updated_dict_unique, updated_dict_abundance)

def filter_after_error_correction(annotations, annkey):
    bad_cdr3 = 0
    num_shm_indels = 0
    bad_primer = 0

    for seq_key in list(annotations[annkey].keys()):
        ann = annotations[annkey][seq_key]
        vprimer = get_vprimer(ann['seq_id']).split("-")[0]
        v_gene =  ann["v_gene"]["fam"]
        if v_gene != vprimer:
            bad_primer += 1
            del annotations[annkey][seq_key]
        elif abstar_cdr3(ann) is None:
            bad_cdr3 += 1
            del annotations[annkey][seq_key]
        elif shm_indels(ann):
            num_shm_indels += 1
            del annotations[annkey][seq_key]

    print(annkey + '_bad_cdr3', bad_cdr3)
    print(annkey + '_num_shm_indels', num_shm_indels)
    print(annkey + '_bad_primer', bad_primer)

def pipeline(f):
    patient = f.split("/")[-1].split(".")[0]
    annotations = load_abstar(f)
    for key in annotations:
        abun = sum([get_abundance(annotations[key][seq_key]['seq_id']) for seq_key in annotations[key]])
        print("initial_" + key, len(annotations[key]), abun)

    run_error_correction(annotations, ['productive'])
    run_error_correction(annotations, ['unproductive'])

    del annotations['stop_in_cdr3']
    del annotations['other']

    filter_after_error_correction(annotations, 'productive')
    filter_after_error_correction(annotations, 'unproductive')
    for key in annotations:
        abun = sum([get_abundance(annotations[key][seq_key]['seq_id']) for seq_key in annotations[key]])
        print('final_', key, len(annotations[key]), abun)

    for key in annotations:
        annotations[key] = list(annotations[key].values())

    save_dir = '/gscratch/stf/zachmon/covid/annotations/abstar/error_corrected/'
    suffix = '_annotations.json'
    #json_save(save_dir + patient + suffix, annotations)
    #print("Saved",save_dir+patient+suffix)

def merge_same_replicate_anns(lin):
    merged_lin = []
    seq_dict = {}
    for ann in lin:
        raw_seq = ann['raw_input'].strip("N")
        time = get_time(ann['seq_id'])
        if (raw_seq,time) not in seq_dict:
            seq_dict[(raw_seq,time)] = []
        seq_dict[(raw_seq,time)].append(ann)
    s_abundance = 0
    for key in seq_dict:
        if len(seq_dict[key]) > 1:
            abundance = sum([get_abundance(ann['seq_id'])
                             for ann in seq_dict[key]])
            s_abundance += abundance
            ann_copy = seq_dict[key][0]
            germ_cdr3 = abstar_naive_cdr3(ann_copy)
            rep_use = get_replicate(ann_copy['seq_id'])
            for a in seq_dict[key][1:]:
                r = get_replicate(a['seq_id'])
                cdr3_a = abstar_naive_cdr3(a)
                if cdr3_a != germ_cdr3:
                    r = get_replicate(a['seq_id'])
                    print("oh no", rep_use, r)
                    print(rep_use,ann_copy['raw_input'],germ_cdr3)
                    print(r, a['raw_input'], cdr3_a)
                    print("\n\n")
            ann_copy['seq_id'] = update_uid(ann_copy['seq_id'], abundance)
            merged_lin.append(ann_copy)
        else:
            merged_lin.append(seq_dict[key][0])
    return merged_lin

def update_uid(uid: str, abundance: int) -> str:
    uid_split = uid.split("|")
    uid_split[3] = "DUPCOUNT=" + str(abundance)
    updated_uid = "|".join(uid_split)
    return updated_uid

def abstar_naive_cdr3(in_a, nt=True, junc=True):
    if nt:
        residue = 'nt'
    else:
        residue = 'aa'
    if junc:
        loc = 'junc_'
    else:
        loc = 'cdr3_'
    try:
        cdr3_start = in_a['vdj_'+residue].index(in_a[loc+residue])
        cdr3_end = cdr3_start + len(in_a[loc+residue])
        return in_a['vdj_germ_'+residue][cdr3_start:cdr3_end]
    except:
        print(in_a['raw_input'].index(in_a['junc_nt']))
        return None

def make_lineages(f):
    patient = f.split("/")[-1].split("_")[0]
    annotations = json_open(f)
    pslc = vjl_slc(annotations['productive'],abstar=True)
    #pslc = merge_replicates(pslc)
    print(len(pslc))
    uslc = vjl_slc(annotations['unproductive'],abstar=True)
    #uslc = merge_replicates(uslc)
    print(len(uslc))
    head_dir = '/gscratch/stf/zachmon/covid/lineages/abstar/error_corrected/'
    suffix = "_lineages_conservative.json"
    save_name = head_dir + patient + suffix
    json_save(save_name, {'productive':pslc, 'unproductive':uslc})
    print("Saved",save_name)

def merge_replicates(lineages):
    out_lins = {}
    for v in lineages:
        for j in lineages[v]:
            for l in lineages[v][j]:
                for cluster_id in lineages[v][j][l]:
                    out_lins[(v,j,l,cluster_id)] = merge_same_replicate_anns(lineages[v][j][l][cluster_id])
    return out_lins
CUTOFF_DICT = {'11': 10, '4': 35,'21': 20,
               '6': 18,'15': 18,'22': 18,
               '9':20,'49':33,'117':10,
               '196':10,'212':10,'271':4,
               '166':10,'187':11,'2':7,
               '42':18,'12':45,'17':25,
               '7':40,'81':0,'82':0,'83':0}

def get_expansion(lineages, patient):
    CUTOFF_DICT = {'11': 10, '4': 35,'21': 20,
                   '6': 18,'15': 18,'22': 18,
                   '9':20,'49':33,'117':10,
                   '196':10,'212':10,'271':4,
                   '166':10,'187':11,'2':7,
                   '42':18,'12':45,'17':25,
                   '7':40,'81':0,'82':0,'83':0}
    df_counts = get_primer_split_counts(lineages, patient, abstar=True,rep=False)['abundance']
    df_fet = fisher_exact_test(df_counts, time_threshold=CUTOFF_DICT[patient], testtype='less')
    exp_lineages = get_expanded_lineages(df_fet, pvalue_thresh=1e-200)
    return exp_lineages

def sort_dict_by_value(x):
    return {k: v
            for k, v in sorted(x.items(),
                               key=lambda item: item[1], reverse=True)}

def get_common_naive_abstar(lineage, nt=True, junc=True):
    naive_cdr3s = []
    for lin_ann in lineage:
        germ_cdr3 = abstar_naive_cdr3(lin_ann, nt=nt, junc=junc)
        naive_cdr3s.append(germ_cdr3)
    counter = sort_dict_by_value(Counter(naive_cdr3s))
    for key in counter:
        if key is None:
            return ""
        if nt:
            if "*" in translate(key):
                continue
            else:
                return key
        else:
            if "*" in key:
                continue
            else:
                return key
    return ""

def create_sonia_input(lin_file):
    in_lineages = json_open(lin_file)['productive']
    patient = lin_file.split("/")[-1].split("_")[0]
    sonia_in = []
    lineages = merge_replicates(in_lineages)
    exp_lineages = get_expansion(lineages, patient)
    for index,li in enumerate(lineages):
        prog_cdr3 = get_common_naive_abstar(lineages[li])
        if prog_cdr3 is None:
            continue
        if len(prog_cdr3) == 0:
            continue
        prog_cdr3_aa = translate(prog_cdr3)
        if prog_cdr3_aa[0] == 'C' and (prog_cdr3_aa[-1] == 'W' or prog_cdr3_aa[-1] == 'F'):
            v = lineages[li][0]["v_gene"]["gene"]
            j = lineages[li][0]["j_gene"]["gene"]
            if li in exp_lineages:
                expanded = True
            else:
                expanded = False
            sonia_in.append((prog_cdr3_aa, v, j, len(lineages[li]), expanded, patient, prog_cdr3))
    df_sonia = pd.DataFrame(sonia_in,columns=['consensus_cdr3',
                                              'v_gene',
                                              'j_gene',
                                              'lineage_size',
                                              'expanded',
                                              'patient',
                                              'nt_cdr3'])
    df_sonia = df_sonia.sort_values(by=['expanded'],ascending=False)
    #  Keep independent nt recombination events only
    df_sonia.drop_duplicates(subset=['v_gene', 'j_gene','consensus_cdr3', 'nt_cdr3'],keep='first')
    sonia_dir = '/gscratch/stf/zachmon/covid/sonia_input/abstar/'
    suffix = '_sonia_input.csv'
    df_sonia.to_csv(sonia_dir + patient + suffix, index=False)
    print("Saved",sonia_dir+patient + suffix)

def main():
    import argparse
    parser = argparse.ArgumentParser(
        description='abstar pipeline')
    parser.add_argument('--in_file', type=str, help='path to abstar file')
    parser.add_argument('--lineages', action='store_true', help='path to abstar file')
    parser.add_argument('--annotations', action='store_true', help='path to abstar file')
    parser.add_argument('--sonia', action='store_true', help='path to abstar file')
    args = parser.parse_args()
    if args.annotations:
       pipeline(args.in_file)
    if args.lineages:
       make_lineages(args.in_file)
    if args.sonia:
        create_sonia_input(args.in_file)

if __name__=='__main__':
    main()

