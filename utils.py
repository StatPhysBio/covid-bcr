import pickle
from os import listdir
from os.path import isfile, join, isdir
import numpy as np
from collections import Counter, OrderedDict
from operator import itemgetter
from jellyfish import hamming_distance

CONST_DATA_DICT = {2:{"samples":['4'],
                       "times":['6'],
                       "severity":"ICU"},
                   6:{"samples":['5','14'],
                      "times":['8','38'],
                      "severity":"Moderate"},
                   11:{"samples":['6','9','15'],
                       "times":['2','15','34'],
                       "severity":"Mild"},
                   15:{"samples":['7','13'],
                       "times":['10','27'],
                       "severity":"Moderate"},
                   21:{"samples":['8','16'],
                       "times":['13','28'],
                       "severity":"Moderate"},
                   22:{"samples":['10','11','18'],
                       "times":['11','16','39'],
                       "severity":"Moderate"},
                   42:{"samples":['17'],
                       "times":['30'],
                       "severity":"ICU"}}

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

def fasta_parse(fasta: str, patient, timepoint, severity, singletons=True):
    seqs = []
    headers = []
    for seq in parse(fasta, 'fasta'):
        uid, cprimer, vprimer, abundance = (x for x in seq.id.split('|'))

        abundance_int = int(abundance.split("=")[-1])
        if not singletons:
            if abundance_int == 1:
                continue

        header_info = [uid, "PATIENT="+str(patient),
                       "TIME=" + str(timepoint), "SEVERITY=" + severity,
                       cprimer, vprimer]
        header = "|".join(header_info)
        header = ">" + header
        headers.append(header)
        seqs.append(str(seq.seq))
    return seqs, headers

def pickle_save(pickle_file, contents):
    with open(pickle_file, 'wb') as handle:
        pickle.dump(contents, handle, protocol=pickle.HIGHEST_PROTOCOL)

def unpickle(pickle_file):
    pickled = pickle.load(open(pickle_file, "rb"))
    return pickled

def get_files(f_dir, suffix):
    files = [f_dir + f
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

def equalize_counters(counter_list):
    all_keys = []
    for c in counter_list:
        all_keys += list(c.keys())
    all_keys = list(set(all_keys))
    for c in counter_list:
        for key in all_keys:
            if key not in c:
                c[key] = 0

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

def get_DUPCOUNT(uid):
    return

def get_dupcounts(uid):
    return int(uid.split("_")[4].replace("dupcounts--",""))

def get_time(uid):
    return int(uid.split("_")[3].replace("time--",""))

def get_primer(uid):
    return uid.split("_")[2].replace("vprimer--","").replace("m", ",")

def get_severity(uid):
    return uid.split("_")[5].replace("severity--","")

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
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }
    protein =""
    if len(seq)%3 == 0:
        for i in range(0, len(seq), 3):
            codon = seq[i:i + 3]
            if "N" in codon:
                print("N in codon!!!!")
                return "("+codon+")"
            protein+= table[codon]
    return protein
