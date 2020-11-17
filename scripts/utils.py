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
import pandas as pd

def trim_bcell_info(df: pd.DataFrame, group):
    """Orients the DataFrame from the csv as a dictionary and deletes duplicate info.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame created from bcell_info.csv.
    group : str
        Specifies which column will be the dictionary keys.

    Returns
    -------
    df_dict : dict
       Dictionary created from DataFrame after removal of duplicate info.
    """

    df_dict = df.drop(group, axis=1).to_dict()

    one_entry = ['organism','age','biomaterial_provider','sex','isolate',
                 'tissue','cell type','disease', 'healthy_state']
    list_entry = ['sample_name','sample_day','replicate']
    dict_entry = ['IgG-nCoV-RBD', 'IgM-nCoV-RBD', 'IgG-SARS-RBD', 'IgG-SARS-NTD',
                  'IgG-nCoV-NTD', 'S2 protein ELISA']

    if group == 'sample_name':
        list_entry.remove(group)
    elif group == 'isolate':
        one_entry.remove(group)

    for key in one_entry:
        df_dict[key] = list(df_dict[key].values())[0]
    for key in dict_entry:
        replacement = {}
        for k,v in df_dict['sample_day'].items():
            replacement[v] = df_dict[key][k]
        df_dict[key] = replacement
    for key in list_entry:
        df_dict[key] = sorted(list(set(df_dict[key].values())))

    df_dict['replicate'] = [int(rep.split(" ")[-1]) for rep in df_dict['replicate']]

    df_dict['severity'] = df_dict['healthy_state']
    del df_dict['healthy_state']

    df_dict['sample day'] = df_dict['sample_day']
    del df_dict['sample_day']

    try:
        df_dict['sample'] = df_dict['sample_name']
        del df_dict['sample_name']
    except:
        pass

    try:
        df_dict['patient'] = df_dict['isolate']
        del df_dict['isolate']
    except:
        pass

    return df_dict

def get_bcell_info(infile: str, group: str='isolate') -> dict:
    """Reads a csv file containg B cell info (following SRA format) and converts it to a dictioanry.

    Parameters
    ----------
    infile : str
        .csv file containing B cell info which is to be read.
    group : str, optional
        Specifies which column will be the dictionary keys.

    Returns
    -------
    df_dict : dict
        Dictionary of B cell info.
    """

    df = pd.read_csv(infile)
    df['isolate'] = df['isolate'].astype(str)
    df_dict = df.groupby(group).apply(lambda dfg: trim_bcell_info(dfg, group)).to_dict()
    return df_dict


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
    """Reads header and sequence information from a FASTA file.

    Parameters
    ----------
    fasta : str
        Path to FASTA.

    Returns
    -------
    sequences : list
        List of sequences in FASTA file.
    headers : list
        List of headers in FASTA file.
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

    dirs = [join(in_dir, d)
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
    """Creates a new FASTA file of the headers and sequences provided.

    Parameters
    ----------
    save_name : str
        Path of FASTA file to be created.
    headers : list
        List of headers of sequences.
    sequences : list
        List of sequences to be saved in FASTA file.

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
