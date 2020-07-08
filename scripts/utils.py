#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""File containing utility functions used by many scripts in this analysis.
    Copyright (C) 2020 Montague, Zachary

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

from collections import Counter, OrderedDict
import csv
import json
from operator import itemgetter
from os import listdir
from os.path import isfile, join, isdir

from Bio.SeqIO import parse

#  Dictionary used when referencing times, severity, plotting binding, etc.
CONST_DATA_DICT = {'H1': {'IgG-nCoV-RBD': [-999.0], 'IgM-nCoV-RBD': [-999.0],
                          'sample day': [0], 'round': 0, 'severity': 'Healthy', 'sample': ['1']},
                   'H2': {'IgG-nCoV-RBD': [-999.0], 'IgM-nCoV-RBD': [-999.0],
                          'sample day': [0], 'round': 0, 'severity': 'Healthy', 'sample': ['2']},
                   'H3': {'IgG-nCoV-RBD': [-999.0], 'IgM-nCoV-RBD': [-999.0],
                          'sample day': [0], 'round': 0, 'severity': 'Healthy', 'sample': ['3']},
                   '19': {'IgG-nCoV-RBD': [0.33], 'IgM-nCoV-RBD': [0.1315],
                          'sample day': [6], 'round': 1, 'severity': 'Severe', 'sample': ['4']},
                   '3': {'IgG-nCoV-RBD': [0.2745, 2.08], 'IgM-nCoV-RBD': [0.361, 1.7305],
                         'sample day': [8, 38], 'round': 1, 'severity': 'Moderate', 'sample': ['5', '14']},
                   '2': {'IgG-nCoV-RBD': [0.799, 1.447, 1.4265], 'IgM-nCoV-RBD': [0.8005, 1.9265, 2.0875],
                         'sample day': [2, 15, 34], 'round': 1, 'severity': 'Mild', 'sample': ['6', '9', '15']},
                   '5': {'IgG-nCoV-RBD': [1.477, 2.028], 'IgM-nCoV-RBD': [1.591, 2.1485],
                         'sample day': [10, 27], 'round': 1, 'severity': 'Moderate', 'sample': ['7', '13']},
                   '6': {'IgG-nCoV-RBD': [1.89, 2.293], 'IgM-nCoV-RBD': [2.904, 2.7875],
                         'sample day': [13, 28], 'round': 1, 'severity': 'Moderate', 'sample': ['8', '16']},
                   '7': {'IgG-nCoV-RBD': [0.1645, 1.218, 2.092], 'IgM-nCoV-RBD': [0.3355, 2.3245, 3.0465],
                         'sample day': [11, 16, 39], 'round': 1, 'severity': 'Moderate', 'sample': ['10', '11', '18']},
                   '18': {'IgG-nCoV-RBD': [2.0025, 1.8605], 'IgM-nCoV-RBD': [1.139, 2.2155],
                          'sample day': [8, 30], 'round': 1, 'severity': 'Severe', 'sample': ['12', '17']},
                   '4': {'IgG-nCoV-RBD': [1.892, 2.276], 'IgM-nCoV-RBD': [1.8885, 2.4285],
                         'sample day': [18, 34], 'round': 2, 'severity': 'Moderate', 'sample': ['19', '23']},
                   '8': {'IgG-nCoV-RBD': [2.3095, 1.9595, 2.023], 'IgM-nCoV-RBD': [2.129, 2.7705, 2.196],
                         'sample day': [14, 32, 37], 'round': 2, 'severity': 'Moderate', 'sample': ['20', '22', '27']},
                   '1': {'IgG-nCoV-RBD': [1.982, 2.0465], 'IgM-nCoV-RBD': [1.823, 1.581],
                         'sample day': [22, 43], 'round': 2, 'severity': 'Mild', 'sample': ['21', '24']},
                   '9': {'IgG-nCoV-RBD': [0.4795, 0.1525, 1.551], 'IgM-nCoV-RBD': [0.904, 0.2125, 2.074],
                         'sample day': [5, 8, 18], 'round': 2, 'severity': 'Moderate', 'sample': ['25', '26', '28']},
                   '16': {'IgG-nCoV-RBD': [1.97, 2.034], 'IgM-nCoV-RBD': [2.107, 1.254],
                          'sample day': [44, 53], 'round': 3, 'severity': 'Severe', 'sample': ['29', '30']},
                   '17': {'IgG-nCoV-RBD': [2.4255, 2.029], 'IgM-nCoV-RBD': [2.848, 2.886],
                          'sample day': [22, 36], 'round': 3, 'severity': 'Severe', 'sample': ['31', '32']},
                   '12': {'IgG-nCoV-RBD': [0.159, 0.152], 'IgM-nCoV-RBD': [1.311, 1.1875],
                          'sample day': [7, 11], 'round': 3, 'severity': 'Moderate', 'sample': ['33', '34']},
                   '13': {'IgG-nCoV-RBD': [0.123, 1.594], 'IgM-nCoV-RBD': [0.1975, 1.0525],
                          'sample day': [7, 13], 'round': 3, 'severity': 'Moderate', 'sample': ['35', '36']},
                   '14': {'IgG-nCoV-RBD': [0.185, 0.308], 'IgM-nCoV-RBD': [0.178, 0.687],
                          'sample day': [3, 7], 'round': 3, 'severity': 'Moderate', 'sample': ['37', '38']},
                   '10': {'IgG-nCoV-RBD': [0.248, 0.7375], 'IgM-nCoV-RBD': [0.3735, 0.726],
                           'sample day': [8, 14], 'round': 3, 'severity': 'Moderate', 'sample': ['39', '40']},
                   '15': {'IgG-nCoV-RBD': [2.0945, 1.9905, 2.0735], 'IgM-nCoV-RBD': [1.0145, 0.8055, 0.6615],
                          'sample day': [39, 41, 71], 'round': 3, 'severity': 'Severe', 'sample': ['41', '42', '43']},
                   '11': {'IgG-nCoV-RBD': [0.1635, 0.518], 'IgM-nCoV-RBD': [0.4795, 1.5075],
                          'sample day': [10, 15], 'round': 3, 'severity': 'Moderate', 'sample': ['44', '45']}}

def csv_to_dict(in_csv_file: str) -> dict:
    """Converts a csv to a dictionary without using pandas.

    Parameters
    ----------
    in_csv_file : str
        Path to csv file.

    Returns
    -------
    csv_dict : dict
        Dictionary of csv data.

    """

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

def fasta_read(fasta: str) -> (list, list):
    """Reads header and sequence information from a fasta file.

    Parameters
    ----------
    fasta : str
        Path to fasta.

    Returns
    -------
    sequences : list
        List of sequences in fasta file.
    headers : list
        List of headers in fasta file.
    """

    sequences = []
    headers = []
    for seq in parse(fasta, 'fasta'):
        headers.append(seq.description)
        sequences.append(str(seq.seq))
    return sequences, headers

def json_save(json_file: str, contents: any) -> None:
    """Saves information to a json file.

    Parameters
    ----------
    json_file : str
        Path to json file to save.
    content : any
        Information to save.

    Returns
    -------
    None
    """

    with open(json_file, 'w') as outfile:
        json.dump(contents, outfile)

def json_open(json_file: str) -> any:
    """Opens a json file and loads the output.

    Parameters
    ----------
    json_file : str
        Path to json file to open.

    Returns
    -------
    contents : any
        Information from the json file.
    """

    with open(json_file) as f:
        contents = json.load(f)
    return contents

def get_files(in_dir: str, keyphrase: str = '') -> list:
    """Returns of list of files with a given keyphrase in a directory.

    Parameters
    ----------
    in_dir : str
        Path to directory containing files of interest.
    keyphrase : str, optional
        String used to filter out other files in a directory.

    Returns
    -------
    files : list
        List of files in a directory.
    """

    files = [join(in_dir, f)
             for f in listdir(in_dir)
             if isfile(join(in_dir, f))
             and ".DS_Store" not in f
             and keyphrase in f]
    files.sort()
    return files

def get_dirs(in_dir: str, keyphrase: str = '') -> list:
    """Returns of list of directories with a given keyphrase in a directory.

    Parameters
    ----------
    in_dir : str
        Path to directory containing directories of interest.
    keyphrase : str, optional
        String used to filter out other directories in a directory.

    Returns
    -------
    dirs : list
        List of dirs in a directory.
    """

    dirs = [join(file_dir, d)
             for d in listdir(in_dir)
             if isdir(join(in_dir, d))
             and ".DS_Store" not in d
             and keyphrase in d]
    dirs.sort()
    return dirs

def create_dir(path: str) -> None:
    """Creates a directory.

    Parameters
    ----------
    path : str
        Path of directory to be created.

    Returns
    -------
    None
    """

    if os.path.exists(path):
        return
    try:
        os.mkdir(path)
    except OSError:
        print ("Creation of the directory %s failed" % path)
    else:
        print ("Successfully created the directory %s " % path)

def create_new_fasta(save_name: str, headers: list, sequences: list) -> None:
    """Creates a new fasta file of the headers and sequences provided.

    Parameters
    ----------
    save_name : str
        Path of fasta file to be created.
    headers : list
        List of headers of sequences.
    sequences : list
        List of sequences to be saved in fasta file.

    Returns
    -------
    None
    """

    with open(save_name, "w") as new_fasta:
        for i, header in enumerate(headers):
            if header[0] != ">":
                header = ">" + header
            new_fasta.write(header + "\n")
            new_fasta.write(sequences[i] + "\n")
    print("Created", save_name)

def sort_dict(unsorted_dict: dict) -> OrderedDict:
    """Sorts dicitionary keys by ascending order.

    Parameters
    ----------
    unsorted_dict : dict
        Dictionary presumed to have its keys not sorted.

    Returns
    -------
    sorted_dict : OrderedDict
        Dictionary with its keys sorted in descending order.
    """

    sorted_dict = OrderedDict(sorted(unsorted_dict.items(), key=itemgetter(0)))
    return sorted_dict

def sort_dict_by_value(unsorted_dict: dict) -> dict:
    """Sorted dictionary keys by value descending order.

    Parameters
    ----------
    unsorted_dict : dict
        Dictionary presumed to have its keys not sorted.

    Returns
    -------
    dict
        Dictionary with its keys sorted according to value descending order.
    """

    return {k: v
            for k, v in sorted(unsorted_dict.items(),
                               key=lambda item: item[1], reverse=True)}


def equalize_counters(counter_list: list, pseudocount: int) -> None:
    """Ensures all counters in the list have the same keys, adding pseudocounts if nonexistent.

    Parameters
    ----------
    counter_list : list
        List of counters.
    pseudocount : int
        Integer to use if key is not present in a counter. Typically use 0 when
        averaging arithmetically and 1 when averaging geometrically.

    Returns
    -------
    None
    """

    all_keys = []
    for c in counter_list:
        all_keys += list(c.keys())
    all_keys = list(set(all_keys))
    for c in counter_list:
        for key in all_keys:
            if key not in c:
                c[key] = pseudocount

def normalize_counter(counter: dict) -> None:
    """Converts counts in counter dictionary to frequencies.

    Parameters
    ----------
    counter : dict
        Dictionary of keys with values being counts.

    Returns
    -------
    None
    """

    normalization = sum(counter.values())
    for key in counter:
        counter[key] /= normalization

def get_patient(header: str) -> str:
    """Obtains the patient from the header.

    Parameters
    ----------
    header : str

    Returns
    -------
    str
        Patient ID.
    """

    return header.split("|")[0].split("=")[-1]

def get_cprimer(header: str) -> str:
    """Obtains the cprimer from the header.

    Parameters
    ----------
    header : str

    Returns
    -------
    str
        cprimer
    """

    return header.split('|')[1].split('=')[-1]

def get_vprimer(header: str) -> str:
    """Obtains the vprimer from the header.

    Parameters
    ----------
    header : str

    Returns
    -------
    str
        vprimer
    """

    return header.split('|')[2].split('=')[-1]

def get_abundance(header: str) -> int:
    """Obtains the abundance counts from the header.

    Parameters
    ----------
    header : str

    Returns
    -------
    int
        abundance counts
    """

    return int(header.split('|')[3].split('=')[-1])

def get_time(header: str) -> int:
    """Obtains the time from the header.

    Parameters
    ----------
    header : str

    Returns
    -------
    int
        time
    """

    return int(header.split('|')[4].split('=')[-1])

def get_severity(header: str) -> str:
    """Obtains the severity from the header.

    Parameters
    ----------
    header : str

    Returns
    -------
    int
        Severity.
    """

    return uid.split("|")[5].split("=")[-1]

def get_replicate(header: str) -> int:
    """Obtains the replicate from the header.

    Parameters
    ----------
    header : str

    Returns
    -------
    int
        replicate
    """

    return int(header.split('|')[6].split('=')[-1])


def remove_N_buffer(sequence: str) -> str:
    """Strips the ambiguous calls (N's) are the beginning and end of a sequence.

    Parameters
    ----------
    sequence : str

    Returns
    -------
    str
        Sequence without any N's at the beginning and end of it.
    """

    return seq.strip("N")
