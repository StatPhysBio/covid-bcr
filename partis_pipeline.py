import pandas as pd
from lineage_fisher_exact_test import *
from error_correct import *
from utils import *
def remove_N_sequences(annotations):
    no_N = []
    for ann in annotations:
        if ann['input_seqs'][0].strip("N").count("N") == 0:
            no_N.append(ann)
    return no_N

def error_correction(annotations):
    ann_dict = {}
    in_h = []
    in_seqs = []
    for ann in annotations:
        ann['unique_ids'][0] = ann['unique_ids'][0]
        key = ann['unique_ids'][0].split("|")[0]
        ann_dict[key] = ann
        in_h.append(ann['unique_ids'][0])
        in_seqs.append(ann['input_seqs'][0].strip("N"))
    #  Error correct
    groups = group_data(in_h,in_seqs)
    for key in groups:
        marg_out_uids, marg_out_seqs, marg_pc = error_correct_marginal(groups[key]['uids'],groups[key]['sequences'])
        tot_out_uids, tot_out_seqs, total_pc = error_correct_total(marg_out_uids,marg_out_seqs,d_tol=2)
        groups[key]['uids'] = tot_out_uids
        groups[key]['sequences'] = tot_out_seqs
    final_uids = []
    final_seqs = []
    for key in groups:
        final_uids += groups[key]['uids']
        final_seqs += groups[key]['sequences']

    #  Get error corrected anns, remove SHM, remove in frame w/ stops
    #  Remove entries whose V gene doesn't match primer
    ec_anns = []
    for out_h in final_uids:
        key = out_h.split("|")[0]
        #  Remove invalid
        if ann_dict[key]['invalid']:
            continue
        #  Remove SHM indels
        if ann_dict[key]['has_shm_indels'][0]:
            continue
        #  Remove in frame with stops
        if ann_dict[key]['in_frames'][0] and ann_dict[key]['stops'][0]:
            continue
        v_gene = ann['v_gene'].split("-")[0]
        vprimer = get_vprimer(ann['unique_ids'][0]).split("-")[0]
        if vprimer == v_gene:
            continue
        ec_anns.append(ann_dict[key])
    return ec_anns

def naive_cdr3(ann):
    cp = ann['codon_positions']
    return ann['naive_seq'][cp['v']:cp['j']+3]

def get_common_naive_partis(lineage):
    naives = [ann['naive_seq'] for ann in lineage]
    cps = [ann['codon_positions'] for ann in lineage]
    naive_cdr3s = [None]*len(naives)
    for i,naive in enumerate(naives):
        naive_cdr3s[i] = naive[cps[i]['v']:cps[i]['j']+3]
    counter = sort_dict_by_value(Counter(naive_cdr3s))
    for key in counter:
        if "*" in translate(key):
            continue
        else:
            return key
    return ""

def remove_wrong_primers(annotations):
    corrected_anns = []
    for ann in annotations:
        v_gene = ann['v_gene'].split("-")[0]
        vprimer = get_vprimer(ann['unique_ids'][0]).split("-")[0]
        if vprimer == v_gene:
            corrected_anns.append(ann)
    return corrected_anns

def update_uid(uid: str, abundance: int) -> str:
    uid_split = uid.split("|")
    uid_split[3] = "DUPCOUNT=" + str(abundance)
    updated_uid = "|".join(uid_split)
    return updated_uid

def merge_same_replicate_anns(lin):
    merged_lin = []
    seq_dict = {}
    for ann in lin:
        raw_seq = ann['input_seqs'][0].strip("N")
        time = get_time(ann['unique_ids'][0])
        if (raw_seq,time) not in seq_dict:
            seq_dict[(raw_seq,time)] = []
        seq_dict[(raw_seq,time)].append(ann)
    s_abundance = 0
    for key in seq_dict:
        if len(seq_dict[key]) > 1:
            abundance = sum([get_abundance(ann['unique_ids'][0])
                             for ann in seq_dict[key]])
            s_abundance += abundance
            ann_copy = seq_dict[key][0]
            germ_cdr3 = naive_cdr3(ann_copy)
            rep_use = get_replicate(ann_copy['unique_ids'][0])
            for a in seq_dict[key][1:]:
                r = get_replicate(a['unique_ids'][0])
                cdr3_a = naive_cdr3(a)
                if cdr3_a != germ_cdr3:
                    r = get_replicate(a['unique_ids'][0])
                    print("oh no", rep_use, r)
                    print(rep_use,ann_copy['input_seqs'][0],germ_cdr3)
                    print(r, a['input_seqs'][0], cdr3_a)
                    print("\n\n")
                else:
                    print("match", rep_use, r)
            ann_copy['unique_ids'][0] = update_uid(ann_copy['unique_ids'][0], abundance)
            merged_lin.append(ann_copy)
        else:
            merged_lin.append(seq_dict[key][0])
    return merged_lin

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
    df_counts = get_primer_split_counts(lineages, patient, rep=False)['abundance']
    df_fet = fisher_exact_test(df_counts, time_threshold=CUTOFF_DICT[patient], testtype='less')
    exp_lineages = get_expanded_lineages(df_fet, pvalue_thresh=1e-200)
    return exp_lineages

def create_sonia_input(patient, lineages, exp_lineages):
    sonia_in = []
    sonia_tuples = []
    for index,li in enumerate(lineages):
        prog_cdr3 = get_common_naive_partis(lineages[li])
        if len(prog_cdr3) == 0:
            continue
        if prog_cdr3[0] == 'C' and (prog_cdr3[-1] == 'W' or prog_cdr3[-1] == 'F'):
            v = lineages[li][0]["v_gene"].split("*")[0]
            j = lineages[li][0]["j_gene"].split("*")[0]
            if li in exp_lineages:
                expanded = True
            else:
                expanded = False
            sonia_in.append((prog_cdr3, v, j, len(lineages[li]), expanded, patient))
    df_sonia = pd.DataFrame(sonia_in,columns=['consensus_cdr3',
                                              'v_gene',
                                              'j_gene',
                                              'lineage_size',
                                              'expanded',
                                              'patient'])
    df_sonia = df_sonia.sort_values(by=['expanded'],ascending=False)
    df_sonia.drop_duplicates(subset=['v_gene', 'j_gene','consensus_cdr3'],keep='first')
    sonia_dir = '/gscratch/stf/zachmon/covid/sonia_input/partis/'
    suffix = '_sonia_input.csv'
    df_sonia.to_csv(sonia_dir + patient + suffix, index=False)
    print("Saved",sonia_dir+patient + suffix)
