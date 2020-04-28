from utils import *
import numpy as np
import operator
import collections
import numpy as np


def initialize_stats_dict(d, label):
    in_frame_no_indels = {"naive seqs": {},
               "input seqs": {},
               "v gene": {},
               "j gene": {},
               "d gene": {},
               "cdr3 length": {},
               "cdr3": {},
               "cdr3 aa":{},
               "cdr3 aa naive":{},
               "vd_mono_freq": {"A":{}, "G":{}, "T":{}, "C":{}},
               "dj_mono_freq": {"A":{}, "G":{}, "T":{}, "C":{}},
               "vd ins": {},
               "dj ins": {},
               "v del": {},
               "j del": {},
               "d5 del": {},
               "d3 del": {},
               "duplicates":{},
               "multiplicities":{},
               "unproductive": []
              }
    in_frame_indels = {"naive seqs": {},
               "input seqs": {},
               "indel reversed seqs": {},
               "indel gapped seqs": {},
               "germline gapped seqs": {},
               "v gene": {},
               "j gene": {},
               "d gene": {},
               "cdr3 length": {},
               "cdr3": {},
               "cdr3 aa":{},
               "cdr3 aa naive":{},
               "vd_mono_freq": {"A":{}, "G":{}, "T":{}, "C":{}},
               "dj_mono_freq": {"A":{}, "G":{}, "T":{}, "C":{}},
               "vd ins": {},
               "dj ins": {},
               "v del": {},
               "j del": {},
               "d5 del": {},
               "d3 del": {},
               "shm ins": {},
               "shm dels": {},
               "shm sum": {},
               "duplicates":{},
               "multiplicities":{},
               "unproductive": []
              }
    out_of_frame_no_indels = {"naive seqs": {},
               "input seqs": {},
               "v gene": {},
               "j gene": {},
               "d gene": {},
               "cdr3 length": {},
               "cdr3": {},
               "cdr3 aa":{},
               "cdr3 aa naive":{},
               "vd_mono_freq": {"A":{}, "G":{}, "T":{}, "C":{}},
               "dj_mono_freq": {"A":{}, "G":{}, "T":{}, "C":{}},
               "vd ins": {},
               "dj ins": {},
               "v del": {},
               "j del": {},
               "d5 del": {},
               "d3 del": {},
               "duplicates":{},
               "multiplicities":{},
               "unproductive": []
              }
    out_of_frame_indels = {"naive seqs": {},
               "input seqs": {},
               "indel reversed seqs": {},
               "indel gapped seqs": {},
               "germline gapped seqs": {},
               "v gene": {},
               "j gene": {},
               "d gene": {},
               "cdr3 length": {},
               "cdr3": {},
               "cdr3 aa":{},
               "cdr3 aa naive":{},
               "vd_mono_freq": {"A":{}, "G":{}, "T":{}, "C":{}},
               "dj_mono_freq": {"A":{}, "G":{}, "T":{}, "C":{}},
               "vd ins": {},
               "dj ins": {},
               "v del": {},
               "j del": {},
               "d5 del": {},
               "d3 del": {},
               "shm ins": {},
               "shm dels": {},
               "shm sum": {},
               "duplicates":{},
               "multiplicities":{},
               "unproductive": []
              }
    d[label] = {"invalid seqs": {},
                "in frame no indels": in_frame_no_indels,
                "in frame indels": in_frame_indels,
                "out of frame no indels": out_of_frame_no_indels,
                "out of frame indels": out_of_frame_indels}

def get_shm_indel_stats(ann):
    num_ins = []
    num_dels = []

    qr_seq_with_gaps=ann['gapped_vdj_nt']
    #qr_gap_positions = find_char(qr_seq_with_gaps, "-")
    #qr_parts = partition_positions(qr_gap_positions)

    gl_seq_with_gaps=ann['gapped_vdj_germ_nt']
    #gl_gap_positions = find_char(gl_seq_with_gaps, "-")
    #gl_parts = partition_positions(gl_gap_positions)

    dels = qr_seq_with_gaps.count("-")
    ins = gl_seq_with_gaps.count("-")
    return dels, ins

def partition_positions(positions):
    partitions = []
    if len(positions) == 1:
        return [positions]
    while positions:
        for i, p in enumerate(positions):
            if p - positions[0] != i:
                partitions.append(positions[:i])
                positions = [x for x in positions
                             if x not in positions[:i]]
                break
            if i == len(positions) - 1:
                partitions.append(positions)
                positions = []
    return partitions

def fill_dict(s_dict, ann, uid):
    s_dict['naive seqs'][uid] = ann['vdj_germ_nt']
    s_dict['input seqs'][uid] = ann['raw_input']

    #  Collect uids with stop codons
    if ann['prod'] == 'no':
        s_dict['unproductive'].append(uid)

    #  Gene statistics
    v_gene = ann['v_gene']['gene']
    j_gene = ann['j_gene']['gene']
    d_gene = ann['d_gene']['gene']
    add_item_to_key(s_dict['v gene'], v_gene, [uid])
    add_item_to_key(s_dict['j gene'], j_gene, [uid])
    add_item_to_key(s_dict['d gene'], d_gene, [uid])

    #  CDR3 statistics
    add_item_to_key(s_dict['cdr3 length'], ann['junc_len'], [uid])
    s_dict['cdr3'][uid] = ann['junc_nt']
    s_dict['cdr3 aa'][uid] = ann['junc_aa']

    #  Mono freq statistics
    vd_ins = ann['junc_nt_breakdown']['n1_nt']
    dj_ins = ann['junc_nt_breakdown']['n2_nt']
    for key in s_dict['vd_mono_freq'].keys():
        s_dict['vd_mono_freq'][key][uid] = vd_ins.count(key)
        s_dict['dj_mono_freq'][key][uid] = dj_ins.count(key)

    #  Insertion statistics
    add_item_to_key(s_dict['vd ins'], len(vd_ins), [uid])
    add_item_to_key(s_dict['dj ins'], len(dj_ins), [uid])

    #  Deletion statistics
    add_item_to_key(s_dict['v del'], ann['exo_trimming']['var_3'], [uid])
    add_item_to_key(s_dict['j del'], ann['exo_trimming']['join_5'], [uid])
    add_item_to_key(s_dict['d5 del'], ann['exo_trimming']['div_5'], [uid])
    add_item_to_key(s_dict['d3 del'], ann['exo_trimming']['div_3'], [uid])


def fill_dict_indels(s_dict, ann, uid):
    #  shm useful seqs
    #  for figuring out other statistics but not needed anymore?
    #s_dict['indel reversed seqs'][uid] = ann['indel_reversed_seqs'][0]
    #s_dict['indel gapped seqs'][uid] = ann['qr_gap_seqs'][0]
    #s_dict['germline gapped seqs'][uid] = ann['gl_gap_seqs'][0]

    #  shm indel statistics
    shm_ins, shm_dels = get_shm_indel_stats(ann)
    add_item_to_key(s_dict['shm ins'], shm_ins, [uid])
    add_item_to_key(s_dict['shm dels'], shm_dels, [uid])
    add_item_to_key(s_dict['shm sum'], shm_ins - shm_dels, [uid])

def get_stats(s_dict, annotations):
    s_dict['total'] = len(list(annotations))
    bad_n = 0
    for j, ann in enumerate(annotations):
        uid = ann['seq_id']
        seq = ann['raw_input']

        #  Throw away sequences with N not put in by partis
        seq_without_N_buffer = remove_N_buffer(seq)
        if seq_without_N_buffer.count("N") > 0:
            bad_n += 1
            continue

        in_frame = ann['junction_in_frame']
        if in_frame:
            frame_key = 'in frame'
        else:
            frame_key = 'out of frame'

        has_shm = (ann['gapped_vdj_nt'] == ann['vdj_nt']) or (ann['gapped_vdj_germ_nt'] == ann['vdj_germ_nt'])
        if has_shm:
            frame_key += ' indels'
        else:
            frame_key += ' no indels'

        fill_dict(s_dict[frame_key], ann, uid)
        if has_shm:
            fill_dict_indels(s_dict[frame_key], ann, uid)

def create_stats_file(in_file):
    file_name = in_file.split("/")[-1]
    patient = file_name.split("_")[0]
    save_dir = in_file.replace(file_name, "")
    save_name = save_dir + patient + "_stats.pickle"
    print("\nBegin unpickling" + in_file)
    full_annotations = unpickle(in_file)
    print("Unpickled")
    s_dict = {}
    initialize_stats_dict(s_dict, patient)
    print("Initialized dict")
    get_stats(s_dict[patient], full_annotations)
    print("Got stats")
    print("Saving stats file", save_name)
    pickle_save(save_name, s_dict)
    print("Pickled stats")

def main():
    import argparse

    parser = argparse.ArgumentParser(
        description='Gets relevant stats from abstar file')
    parser.add_argument('--infile', type=str, help='path to abstar pickled file')
    args = parser.parse_args()
    create_stats_file(args.infile)

if __name__ == '__main__':
    main()
