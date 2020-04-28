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
               "stops": []
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
               "stops": []
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
               "stops": []
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
               "stops": []
              }
    d[label] = {"invalid seqs": {},
                "in frame no indels": in_frame_no_indels,
                "in frame indels": in_frame_indels,
                "out of frame no indels": out_of_frame_no_indels,
                "out of frame indels": out_of_frame_indels}

def get_shm_indel_stats(ann):
    num_ins = []
    num_dels = []

    qr_seq_with_gaps=ann['qr_gap_seqs'][0]
    qr_gap_positions = find_char(qr_seq_with_gaps, ".")
    #qr_parts = partition_positions(qr_gap_positions)

    gl_seq_with_gaps=ann['gl_gap_seqs'][0]
    gl_gap_positions = find_char(gl_seq_with_gaps, ".")
    #gl_parts = partition_positions(gl_gap_positions)

    dels = len(qr_gap_positions)
    ins = len(gl_gap_positions)
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
    s_dict['naive seqs'][uid] = ann['naive_seq']
    s_dict['input seqs'][uid] = ann['input_seqs'][0]

    #  Collect uids with stop codons
    if ann['stops'][0]:
        s_dict['stops'].append(uid)

    #  Gene statistics
    v_gene = ann['v_gene'].split("*")[0]
    j_gene = ann['j_gene'].split("*")[0]
    d_gene = ann['d_gene'].split("*")[0]
    add_item_to_key(s_dict['v gene'], v_gene, [uid])
    add_item_to_key(s_dict['j gene'], j_gene, [uid])
    add_item_to_key(s_dict['d gene'], d_gene, [uid])

    #  CDR3 statistics
    add_item_to_key(s_dict['cdr3 length'], ann['cdr3_length'], [uid])
    s_dict['cdr3'][uid] = ann['cdr3_seqs'][0]
    if ann['cdr3_seqs'][0].count("N") == 0:
        s_dict['cdr3 aa'][uid] = translate(ann['cdr3_seqs'][0])
        v_gene_in = ann['codon_positions']['v']
        cdr3_len = ann['cdr3_length']
        s_dict['cdr3 aa naive'][uid] = ann['naive_seq_aa'][int(v_gene_in/3):int((v_gene_in+cdr3_len)/3)]

    if ann['cdr3_length'] != len(ann['cdr3_seqs'][0]):
        print(ann)
    #  Mono freq statistics
    for key in s_dict['vd_mono_freq'].keys():
        s_dict['vd_mono_freq'][key][uid] = ann['vd_insertion'].count(key)
        s_dict['dj_mono_freq'][key][uid] = ann['dj_insertion'].count(key)

    #  Insertion statistics
    add_item_to_key(s_dict['vd ins'], len(ann['vd_insertion']), [uid])
    add_item_to_key(s_dict['dj ins'], len(ann['dj_insertion']), [uid])

    #  Deletion statistics
    add_item_to_key(s_dict['v del'], ann['v_3p_del'], [uid])
    add_item_to_key(s_dict['j del'], ann['j_5p_del'], [uid])
    add_item_to_key(s_dict['d5 del'], ann['d_5p_del'], [uid])
    add_item_to_key(s_dict['d3 del'], ann['d_3p_del'], [uid])


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
        uid = ann['unique_ids'][0]
        seq = ann['input_seqs'][0]

        #  Nothing in this entry
        if ann['invalid'] == True:
            s_dict['invalid seqs'][uid] = seq
            continue

        #  Throw away sequences with N not put in by partis
        seq_without_N_buffer = remove_N_buffer(seq)
        if seq_without_N_buffer.count("N") > 0:
            bad_n += 1
            continue

        in_frame = ann['in_frames'][0]
        if in_frame:
            frame_key = 'in frame'
        else:
            frame_key = 'out of frame'

        has_shm = ann['has_shm_indels'][0]
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
    save_dir = in_file.replace(filename, "")
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
        description='Gets relevant stats from partis file')
    parser.add_argument('--infile', type=str, help='path to partis pickled file')
    args = parser.parse_args()
    create_stats_file(args.infile)

if __name__ == '__main__':
    main()


