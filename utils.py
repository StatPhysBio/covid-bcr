import pickle
from os import listdir
from os.path import isfile, join, isdir
import numpy as np
from collections import Counter, OrderedDict
from operator import itemgetter
from jellyfish import hamming_distance
import csv
from Bio.SeqIO import parse
import json

CONST_SAMPLE_DICT = {'1': {'IgG-nCoV-RBD': -999.0, 'IgM-nCoV-RBD': -999.0, 'sample day': 0, 'round': 0, 'severity': 'Healthy', 'patient': '81'},
                     '2': {'IgG-nCoV-RBD': -999.0, 'IgM-nCoV-RBD': -999.0, 'sample day': 0, 'round': 0, 'severity': 'Healthy', 'patient': '82'},
                     '3': {'IgG-nCoV-RBD': -999.0, 'IgM-nCoV-RBD': -999.0, 'sample day': 0, 'round': 0, 'severity': 'Healthy', 'patient': '83'},
                     '4': {'IgG-nCoV-RBD': 0.33, 'IgM-nCoV-RBD': 0.1315, 'sample day': 6, 'round': 1, 'severity': 'Severe', 'patient': '2'},
                     '5': {'IgG-nCoV-RBD': 0.2745, 'IgM-nCoV-RBD': 0.361, 'sample day': 8, 'round': 1, 'severity': 'Moderate', 'patient': '6'},
                     '6': {'IgG-nCoV-RBD': 0.799, 'IgM-nCoV-RBD': 0.8005, 'sample day': 2, 'round': 1, 'severity': 'Mild', 'patient': '11'},
                     '7': {'IgG-nCoV-RBD': 1.477, 'IgM-nCoV-RBD': 1.591, 'sample day': 10, 'round': 1, 'severity': 'Moderate', 'patient': '15'},
                     '8': {'IgG-nCoV-RBD': 1.89, 'IgM-nCoV-RBD': 2.904, 'sample day': 13, 'round': 1, 'severity': 'Moderate', 'patient': '21'},
                     '9': {'IgG-nCoV-RBD': 1.447, 'IgM-nCoV-RBD': 1.9265, 'sample day': 15, 'round': 1, 'severity': 'Mild', 'patient': '11'},
                     '10': {'IgG-nCoV-RBD': 0.1645, 'IgM-nCoV-RBD': 0.3355, 'sample day': 11, 'round': 1, 'severity': 'Moderate', 'patient': '22'},
                     '11': {'IgG-nCoV-RBD': 1.218, 'IgM-nCoV-RBD': 2.3245, 'sample day': 16, 'round': 1, 'severity': 'Moderate', 'patient': '22'},
                     '12': {'IgG-nCoV-RBD': 2.0025, 'IgM-nCoV-RBD': 1.139, 'sample day': 8, 'round': 1, 'severity': 'Severe', 'patient': '42'},
                     '13': {'IgG-nCoV-RBD': 2.028, 'IgM-nCoV-RBD': 2.1485, 'sample day': 27, 'round': 1, 'severity': 'Moderate', 'patient': '15'},
                     '14': {'IgG-nCoV-RBD': 2.08, 'IgM-nCoV-RBD': 1.7305, 'sample day': 38, 'round': 1, 'severity': 'Moderate', 'patient': '6'},
                     '15': {'IgG-nCoV-RBD': 1.4265, 'IgM-nCoV-RBD': 2.0875, 'sample day': 34, 'round': 1, 'severity': 'Mild', 'patient': '11'},
                     '16': {'IgG-nCoV-RBD': 2.293, 'IgM-nCoV-RBD': 2.7875, 'sample day': 28, 'round': 1, 'severity': 'Moderate', 'patient': '21'},
                     '17': {'IgG-nCoV-RBD': 1.8605, 'IgM-nCoV-RBD': 2.2155, 'sample day': 30, 'round': 1, 'severity': 'Severe', 'patient': '42'},
                     '18': {'IgG-nCoV-RBD': 2.092, 'IgM-nCoV-RBD': 3.0465, 'sample day': 39, 'round': 1, 'severity': 'Moderate', 'patient': '22'},
                     '19': {'IgG-nCoV-RBD': 1.892, 'IgM-nCoV-RBD': 1.8885, 'sample day': 18, 'round': 2, 'severity': 'Moderate', 'patient': '9'},
                     '20': {'IgG-nCoV-RBD': 2.3095, 'IgM-nCoV-RBD': 2.129, 'sample day': 14, 'round': 2, 'severity': 'Moderate', 'patient': '49'},
                     '21': {'IgG-nCoV-RBD': 1.982, 'IgM-nCoV-RBD': 1.823, 'sample day': 22, 'round': 2, 'severity': 'Mild', 'patient': '4'},
                     '22': {'IgG-nCoV-RBD': 1.9595, 'IgM-nCoV-RBD': 2.7705, 'sample day': 32, 'round': 2, 'severity': 'Moderate', 'patient': '49'},
                     '23': {'IgG-nCoV-RBD': 2.276, 'IgM-nCoV-RBD': 2.4285, 'sample day': 34, 'round': 2, 'severity': 'Moderate', 'patient': '9'},
                     '24': {'IgG-nCoV-RBD': 2.0465, 'IgM-nCoV-RBD': 1.581, 'sample day': 43, 'round': 2, 'severity': 'Mild', 'patient': '4'},
                     '25': {'IgG-nCoV-RBD': 0.4795, 'IgM-nCoV-RBD': 0.904, 'sample day': 5, 'round': 2, 'severity': 'Moderate', 'patient': '117'},
                     '26': {'IgG-nCoV-RBD': 0.1525, 'IgM-nCoV-RBD': 0.2125, 'sample day': 8, 'round': 2, 'severity': 'Moderate', 'patient': '117'},
                     '27': {'IgG-nCoV-RBD': 2.023, 'IgM-nCoV-RBD': 2.196, 'sample day': 37, 'round': 2, 'severity': 'Moderate', 'patient': '49'},
                     '28': {'IgG-nCoV-RBD': 1.551, 'IgM-nCoV-RBD': 2.074, 'sample day': 18, 'round': 2, 'severity': 'Moderate', 'patient': '117'},
                     '29': {'IgG-nCoV-RBD': 1.97, 'IgM-nCoV-RBD': 2.107, 'sample day': 44, 'round': 3, 'severity': 'Severe', 'patient': '12'},
                     '30': {'IgG-nCoV-RBD': 2.034, 'IgM-nCoV-RBD': 1.254, 'sample day': 53, 'round': 3, 'severity': 'Severe', 'patient': '12'},
                     '31': {'IgG-nCoV-RBD': 2.4255, 'IgM-nCoV-RBD': 2.848, 'sample day': 22, 'round': 3, 'severity': 'Severe', 'patient': '17'},
                     '32': {'IgG-nCoV-RBD': 2.029, 'IgM-nCoV-RBD': 2.886, 'sample day': 36, 'round': 3, 'severity': 'Severe', 'patient': '17'},
                     '33': {'IgG-nCoV-RBD': 0.159, 'IgM-nCoV-RBD': 1.311, 'sample day': 7, 'round': 3, 'severity': 'Moderate', 'patient': '196'},
                     '34': {'IgG-nCoV-RBD': 0.152, 'IgM-nCoV-RBD': 1.1875, 'sample day': 11, 'round': 3, 'severity': 'Moderate', 'patient': '196'},
                     '35': {'IgG-nCoV-RBD': 0.123, 'IgM-nCoV-RBD': 0.1975, 'sample day': 7, 'round': 3, 'severity': 'Moderate', 'patient': '212'},
                     '36': {'IgG-nCoV-RBD': 1.594, 'IgM-nCoV-RBD': 1.0525, 'sample day': 13, 'round': 3, 'severity': 'Moderate', 'patient': '212'},
                     '37': {'IgG-nCoV-RBD': 0.185, 'IgM-nCoV-RBD': 0.178, 'sample day': 3, 'round': 3, 'severity': 'Moderate', 'patient': '271'},
                     '38': {'IgG-nCoV-RBD': 0.308, 'IgM-nCoV-RBD': 0.687, 'sample day': 7, 'round': 3, 'severity': 'Moderate', 'patient': '271'},
                     '39': {'IgG-nCoV-RBD': 0.248, 'IgM-nCoV-RBD': 0.3735, 'sample day': 8, 'round': 3, 'severity': 'Moderate', 'patient': '166'},
                     '40': {'IgG-nCoV-RBD': 0.7375, 'IgM-nCoV-RBD': 0.726, 'sample day': 14, 'round': 3, 'severity': 'Moderate', 'patient': '166'},
                     '41': {'IgG-nCoV-RBD': 2.0945, 'IgM-nCoV-RBD': 1.0145, 'sample day': 39, 'round': 3, 'severity': 'Severe', 'patient': '7'},
                     '42': {'IgG-nCoV-RBD': 1.9905, 'IgM-nCoV-RBD': 0.8055, 'sample day': 41, 'round': 3, 'severity': 'Severe', 'patient': '7'},
                     '43': {'IgG-nCoV-RBD': 2.0735, 'IgM-nCoV-RBD': 0.6615, 'sample day': 71, 'round': 3, 'severity': 'Severe', 'patient': '7'},
                     '44': {'IgG-nCoV-RBD': 0.1635, 'IgM-nCoV-RBD': 0.4795, 'sample day': 10, 'round': 3, 'severity': 'Moderate', 'patient': '187'},
                     '45': {'IgG-nCoV-RBD': 0.518, 'IgM-nCoV-RBD': 1.5075, 'sample day': 15, 'round': 3, 'severity': 'Moderate', 'patient': '187'}}
CONST_DATA_DICT = {'81': {'IgG-nCoV-RBD': [-999.0], 'IgM-nCoV-RBD': [-999.0], 'sample day': [0], 'round': 0, 'severity': 'Healthy', 'sample': ['1']},
                   '82': {'IgG-nCoV-RBD': [-999.0], 'IgM-nCoV-RBD': [-999.0], 'sample day': [0], 'round': 0, 'severity': 'Healthy', 'sample': ['2']},
                   '83': {'IgG-nCoV-RBD': [-999.0], 'IgM-nCoV-RBD': [-999.0], 'sample day': [0], 'round': 0, 'severity': 'Healthy', 'sample': ['3']},
                   '2': {'IgG-nCoV-RBD': [0.33], 'IgM-nCoV-RBD': [0.1315], 'sample day': [6], 'round': 1, 'severity': 'Severe', 'sample': ['4']},
                   '6': {'IgG-nCoV-RBD': [0.2745, 2.08], 'IgM-nCoV-RBD': [0.361, 1.7305], 'sample day': [8, 38], 'round': 1, 'severity': 'Moderate', 'sample': ['5', '14']},
                   '11': {'IgG-nCoV-RBD': [0.799, 1.447, 1.4265], 'IgM-nCoV-RBD': [0.8005, 1.9265, 2.0875], 'sample day': [2, 15, 34], 'round': 1, 'severity': 'Mild', 'sample': ['6', '9', '15']},
                   '15': {'IgG-nCoV-RBD': [1.477, 2.028], 'IgM-nCoV-RBD': [1.591, 2.1485], 'sample day': [10, 27], 'round': 1, 'severity': 'Moderate', 'sample': ['7', '13']},
                   '21': {'IgG-nCoV-RBD': [1.89, 2.293], 'IgM-nCoV-RBD': [2.904, 2.7875], 'sample day': [13, 28], 'round': 1, 'severity': 'Moderate', 'sample': ['8', '16']},
                   '22': {'IgG-nCoV-RBD': [0.1645, 1.218, 2.092], 'IgM-nCoV-RBD': [0.3355, 2.3245, 3.0465], 'sample day': [11, 16, 39], 'round': 1, 'severity': 'Moderate', 'sample': ['10', '11', '18']},
                   '42': {'IgG-nCoV-RBD': [2.0025, 1.8605], 'IgM-nCoV-RBD': [1.139, 2.2155], 'sample day': [8, 30], 'round': 1, 'severity': 'Severe', 'sample': ['12', '17']},
                   '9': {'IgG-nCoV-RBD': [1.892, 2.276], 'IgM-nCoV-RBD': [1.8885, 2.4285], 'sample day': [18, 34], 'round': 2, 'severity': 'Moderate', 'sample': ['19', '23']},
                   '49': {'IgG-nCoV-RBD': [2.3095, 1.9595, 2.023], 'IgM-nCoV-RBD': [2.129, 2.7705, 2.196], 'sample day': [14, 32, 37], 'round': 2, 'severity': 'Moderate', 'sample': ['20', '22', '27']},
                   '4': {'IgG-nCoV-RBD': [1.982, 2.0465], 'IgM-nCoV-RBD': [1.823, 1.581], 'sample day': [22, 43], 'round': 2, 'severity': 'Mild', 'sample': ['21', '24']},
                   '117': {'IgG-nCoV-RBD': [0.4795, 0.1525, 1.551], 'IgM-nCoV-RBD': [0.904, 0.2125, 2.074], 'sample day': [5, 8, 18], 'round': 2, 'severity': 'Moderate', 'sample': ['25', '26', '28']},
                   '12': {'IgG-nCoV-RBD': [1.97, 2.034], 'IgM-nCoV-RBD': [2.107, 1.254], 'sample day': [44, 53], 'round': 3, 'severity': 'Severe', 'sample': ['29', '30']},
                   '17': {'IgG-nCoV-RBD': [2.4255, 2.029], 'IgM-nCoV-RBD': [2.848, 2.886], 'sample day': [22, 36], 'round': 3, 'severity': 'Severe', 'sample': ['31', '32']},
                   '196': {'IgG-nCoV-RBD': [0.159, 0.152], 'IgM-nCoV-RBD': [1.311, 1.1875], 'sample day': [7, 11], 'round': 3, 'severity': 'Moderate', 'sample': ['33', '34']},
                   '212': {'IgG-nCoV-RBD': [0.123, 1.594], 'IgM-nCoV-RBD': [0.1975, 1.0525], 'sample day': [7, 13], 'round': 3, 'severity': 'Moderate', 'sample': ['35', '36']},
                   '271': {'IgG-nCoV-RBD': [0.185, 0.308], 'IgM-nCoV-RBD': [0.178, 0.687], 'sample day': [3, 7], 'round': 3, 'severity': 'Moderate', 'sample': ['37', '38']},
                   '166': {'IgG-nCoV-RBD': [0.248, 0.7375], 'IgM-nCoV-RBD': [0.3735, 0.726], 'sample day': [8, 14], 'round': 3, 'severity': 'Moderate', 'sample': ['39', '40']},
                   '7': {'IgG-nCoV-RBD': [2.0945, 1.9905, 2.0735], 'IgM-nCoV-RBD': [1.0145, 0.8055, 0.6615], 'sample day': [39, 41, 71], 'round': 3, 'severity': 'Severe', 'sample': ['41', '42', '43']},
                   '187': {'IgG-nCoV-RBD': [0.1635, 0.518], 'IgM-nCoV-RBD': [0.4795, 1.5075], 'sample day': [10, 15], 'round': 3, 'severity': 'Moderate', 'sample': ['44', '45']}}

partis_v_genes = ['IGHV1-18', 'IGHV1-2', 'IGHV1-24', 'IGHV1-3',
                  'IGHV1-46', 'IGHV1-58', 'IGHV1-69', 'IGHV1-8',
                  'IGHV2-26', 'IGHV2-5', 'IGHV2-70', 'IGHV3-11',
                  'IGHV3-13', 'IGHV3-15', 'IGHV3-20', 'IGHV3-21',
                  'IGHV3-23', 'IGHV3-30', 'IGHV3-30-3', 'IGHV3-33',
                  'IGHV3-43', 'IGHV3-48', 'IGHV3-49', 'IGHV3-53',
                  'IGHV3-64', 'IGHV3-66', 'IGHV3-7', 'IGHV3-72',
                  'IGHV3-73', 'IGHV3-74', 'IGHV3-9', 'IGHV4-30-2',
                  'IGHV4-30-4', 'IGHV4-31', 'IGHV4-34', 'IGHV4-38-2',
                  'IGHV4-39', 'IGHV4-4', 'IGHV4-59', 'IGHV4-61',
                  'IGHV5-10-1', 'IGHV5-51', 'IGHV6-1', 'IGHV7-4-1']

partis_j_genes = ['IGHJ1', 'IGHJ2', 'IGHJ3', 'IGHJ4', 'IGHJ5', 'IGHJ6']

abstar_v_genes = ['IGHV1-18', 'IGHV1-2', 'IGHV1-24', 'IGHV1-3',
                  'IGHV1-45', 'IGHV1-46', 'IGHV1-58', 'IGHV1-69',
                  'IGHV1-69-2', 'IGHV1-8', 'IGHV2-26', 'IGHV2-5',
                  'IGHV2-70', 'IGHV2-70D', 'IGHV3-11', 'IGHV3-13',
                  'IGHV3-15', 'IGHV3-20', 'IGHV3-21', 'IGHV3-23',
                  'IGHV3-30', 'IGHV3-30-3', 'IGHV3-33', 'IGHV3-43',
                  'IGHV3-43D', 'IGHV3-48', 'IGHV3-49', 'IGHV3-53',
                  'IGHV3-64', 'IGHV3-64D', 'IGHV3-66', 'IGHV3-7',
                  'IGHV3-72', 'IGHV3-73', 'IGHV3-74', 'IGHV3-9',
                  'IGHV3-NL1', 'IGHV4-28', 'IGHV4-30-2', 'IGHV4-30-4',
                  'IGHV4-31', 'IGHV4-34', 'IGHV4-38-2', 'IGHV4-39',
                  'IGHV4-4', 'IGHV4-59', 'IGHV4-61', 'IGHV5-10-1',
                  'IGHV5-51', 'IGHV6-1', 'IGHV7-4-1']

abstar_j_genes = ['IGHJ1', 'IGHJ2', 'IGHJ3', 'IGHJ4', 'IGHJ5', 'IGHJ6']

def csv_to_dict(in_csv_file):
    csv_dict = {}

    with open(in_csv_file, newline='',encoding='utf-8-sig') as csvfile:
        reader = csv.DictReader(csvfile)
        for index,ordered_dict in enumerate(reader):
            if list(ordered_dict.values())==['']*len(ordered_dict):
                break
            for key in ordered_dict.keys():
                if index == 0:
                    csv_dict[key]=[]
                csv_dict[key].append(ordered_dict[key])
    return csv_dict

def fasta_read(fasta: str):
    seqs = []
    headers = []
    for seq in parse(fasta, 'fasta'):
        headers.append(seq.id)
        seqs.append(str(seq.seq))
    return seqs, headers

def fasta_parse(fasta: str, patient='', timepoint='',
                severity='', replicate='',
                oneline_collapsed=False, singletons=True):
    seqs = []
    headers = []

    if oneline_collapsed:
        filename = fasta.split("/")[-1]
        sample = filename.split("_")[0].replace("S","")
        if "-" in sample:
            replicate = sample.split("-")[-1]
        else:
            replicate = '1'

    for seq in parse(fasta, 'fasta'):
        uid, cprimer, vprimer, abundance = (x for x in seq.id.split('|'))

        abundance_int = int(abundance.split("=")[-1])
        if not singletons:
            if abundance_int == 1:
                continue

        if oneline_collapsed:
            sample = sample.split("-")[0]
            header_info = [uid+"="+CONST_SAMPLE_DICT[sample]['patient'],
                           cprimer, vprimer,abundance,
                           "TIME=" + str(CONST_SAMPLE_DICT[sample]['sample day']),
                           "SEVERITY=" + str(CONST_SAMPLE_DICT[sample]['severity']),
                           "REPLICATE=" + replicate]
        else:
            header_info = [uid+"="+patient,
                           cprimer, vprimer,abundance,
                           "TIME=" + timepoint,
                           "SEVERITY=" + severity,
                           "REPLICATE=" + replicate]

        header = "|".join(header_info)
        header = ">" + header
        headers.append(header)
        seqs.append(str(seq.seq))
    return seqs, headers

def json_save(json_file, contents):
    with open(json_file, 'w') as outfile:
        json.dump(contents, outfile)

def json_open(json_file):
    with open(json_file) as f:
        contents = json.load(f)
    return contents

def pickle_save(pickle_file, contents):
    with open(pickle_file, 'wb') as handle:
        pickle.dump(contents, handle, protocol=pickle.HIGHEST_PROTOCOL)

def unpickle(pickle_file):
    pickled = pickle.load(open(pickle_file, "rb"))
    return pickled

def get_files(f_dir, suffix):
    files = [join(f_dir, f)
             for f in listdir(f_dir)
             if isfile(join(f_dir, f))
             and ".DS_Store" not in f
             and suffix in f]
    files.sort()
    return files

def create_dir(path):
    if os.path.exists(path):
        return
    try:
        os.mkdir(path)
    except OSError:
        print ("Creation of the directory %s failed" % path)
    else:
        print ("Successfully created the directory %s " % path)

def create_new_fasta(save_name, headers, sequences):
    with open(save_name, "w") as new_fasta:
        for i, header in enumerate(headers):
            new_fasta.write(header + "\n")
            new_fasta.write(sequences[i] + "\n")

def add_item_to_key(d, key, item):
    if key in d:
        d[key] += item
    else:
        d[key] = item

def update_dictionary(d, temp_d):
    for key in temp_d:
        if key in d:
            d[key] += temp_d[key]
        else:
            d[key] = temp_d[key]

def sort_dict(unsorted_dict):
    D = OrderedDict(sorted(unsorted_dict.items(), key=itemgetter(0)))
    return D

def get_data_counter(data, multiplicity=False):
    data_counter = {}
    if not multiplicity:
        for key in data:
            data_counter[key] = len(data[key])
        return data_counter

    for key in data:
        data_counter[key] = 0
        for uid in data[key]:
            data_counter[key] += get_dupcounts(uid)
    return data_counter

def equalize_counters(counter_list, pseudocount):
    all_keys = []
    for c in counter_list:
        all_keys += list(c.keys())
    all_keys = list(set(all_keys))
    for c in counter_list:
        for key in all_keys:
            if key not in c:
                c[key] = pseudocount

def normalize_counter(data_counter):
    normalization = get_total_counts(data_counter)
    for key in data_counter:
        data_counter[key] /= normalization
    return data_counter

def get_total_counts(data_counter):
    return sum(data_counter.values())

def split_dict_by_time(in_dict, time_threshold=None):
    time_dicts = {}
    time_el_dict = {"data": {"early":{},
                             "late": {}},
                    "times":{"early": [],
                             "late": []}}
    for key in in_dict:
        for item in in_dict[key]:
            time = get_time(item)
            if time not in time_dicts:
                time_dicts[time] = {}
            if key not in time_dicts[time]:
                time_dicts[time][key] = []
            time_dicts[time][key].append(item)
    if time_threshold != None:
        for time in time_dicts:
            if time < time_threshold:
                time_el_dict['times']['early'].append(time)
                update_dictionary(time_el_dict['data']['early'], time_dicts[time])
            else:
                time_el_dict['times']['late'].append(time)
                update_dictionary(time_el_dict['data']['late'], time_dicts[time])
        return time_el_dict
    return time_dicts

def split_dict_by_primer(in_dict):
    primer_dict = {}
    for key in in_dict:
        for uid in in_dict[key]:
            primer = get_primer(uid)
            if primer not in primer_dict:
                primer_dict[primer] = {}
            if key not in primer_dict[primer]:
                primer_dict[primer][key] = []
            primer_dict[primer][key].append(uid)
    sorted_primer_dict = sort_dict(primer_dict)
    primer_list = list(sorted_primer_dict.keys())
    bad_keys = [p for p in primer_list if "," in p]
    for p in bad_keys:
        del sorted_primer_dict[p]
    return sorted_primer_dict

def convert_header(header):
    header_split = header.split("_")
    seq_id = header_split[0]
    patient_num = header_split[1].split("--")[-1]
    vprimer = header_split[2].split("--")[-1]
    time = header_split[3].split("--")[-1]
    abundance = header_split[4].split("--")[-1]
    severity = header_split[5].split("--")[-1]
    return "|".join([seq_id + "-" + patient_num,
                     'CPRIMER=CHG-R',
                     'VPRIMER='+vprimer,
                     'DUPCOUNT='+abundance,
                     'TIME='+time,
                     'SEVERITY='+severity])

def get_patient(uid: str) -> int:
    return int(uid.split("|")[0].split("=")[-1])

def get_cprimer(uid: str) -> str:
    return uid.split("|")[1].split("=")[-1]

def get_vprimer(uid: str) -> str:
    return uid.split("|")[2].split("=")[-1]

def get_abundance(uid: str) -> int:
    return int(uid.split("|")[3].split("=")[-1])

def get_time(uid: str) -> int:
    return int(uid.split("|")[4].split("=")[-1])

def get_severity(uid: str) -> str:
    return uid.split("|")[5].split("=")[-1]

def get_replicate(uid: str) -> str:
    return uid.split("|")[6].split("=")[-1]

def find_char(string, c):
    return [pos for pos, char in enumerate(string) if char == c]

def get_n_at_end(n_list, s_len):
    last_index = s_len - 1
    n_list_rev = list(reversed(n_list))
    end_buffer = []
    n_list_rev_copy = n_list_rev.copy()
    for i, loc in enumerate(n_list_rev):
        if loc == last_index - i:
            end_buffer.append(last_index - i)
    return list(reversed(end_buffer))

def get_n_at_beginning(n_list):
    beginning_buffer = []
    for i, loc in enumerate(n_list):
        if loc == i:
            beginning_buffer.append(i)
    return beginning_buffer

def get_num_n_at_start_and_end(seq):
    n_locs = find_char(seq, "N")
    num_start = len(get_n_at_beginning(n_locs))
    num_end = len(get_n_at_end(n_locs, len(seq)))
    return num_start, num_end

def remove_N_buffer(seq):
    n_locs = find_char(seq, "N")
    beg_n = get_n_at_beginning(n_locs)
    n_locs = [n for n in n_locs if n not in beg_n]
    end_n = get_n_at_end(n_locs, len(seq))
    lower_bound = 0
    if beg_n:
        lower_bound = np.max(beg_n) +1
    upper_bound = len(seq)
    if end_n:
        upper_bound = np.min(end_n)
    return seq[lower_bound:upper_bound]

def fix_N_padding(seq_to_fix, new_num_n_start, new_num_n_end):
    seq_no_pad = remove_N_buffer(seq_to_fix)
    fixed_seq = "".join("N"*new_num_n_start)+seq_no_pad + "".join("N"*new_num_n_end)
    return fixed_seq

def translate(seq):

    table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
        'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
    }
    protein =""
    if len(seq)%3 == 0:
        for i in range(0, len(seq), 3):
            codon = seq[i:i + 3]
            if "N" in codon:
                protein += 'X'
            else:
                protein+= table[codon]
    return protein
