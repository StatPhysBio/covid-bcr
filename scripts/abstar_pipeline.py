#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Script to run the entire pipeline performed on sequences used in paper.
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

import pandas as pd

from error_correct import error_correct_marginal, error_correct_total, group_data
from utils import *
from vjl_slc import vjl_slc

def translate(sequence:str) -> str:
    """Takes a nucleotide sequence and translates it into a protein.

    Parameters
    ----------
    sequence : str
        String of nucleotide residues.

    Returns
    -------
    protein : str
        String of amino acids.
    """

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
    if len(sequence) % 3 == 0:
        for i in range(0, len(sequence), 3):
            codon = sequence[i:i + 3]
            if "N" in codon:
                protein += 'X'
            else:
                protein += table[codon]
    return protein


def load_abstar(abstar_file: str) -> dict:
    """Loads and performs preliminary filtering on abstar output.

    Some json dictionaries from abstar output may be corrupted and therefore
    difficult to load. The abstar file is read dictionary-by-dictionary
    instead of all at once. Annotations which are not heavy chain,
    missing CDR3 when aligned, or contain an N in the input sequence
    are skipped. All passed sequences are saved to a dictionary
    which splits them into productive, unproductive (out-of-frame exclusive),
    stop codon in CDR3, and all other cases. If one were to perform error
    correction more liberally, including the annotations in "other" is
    more permissble than using annotations which show that there's a
    stop codon in the CDR3.

    Parameters
    ----------
    abstar_file : str
        Path to abstar .json output.

    Returns
    -------
    annotations : dict
        Dictionary of directiony of annotations split by types of productivity.
        keys of the first dictionary are ['productive', 'unproductive',
        'stop_in_cdr3', 'other']. Keys in the second dictionary come
        from the sequence header/id for each annotation. This makes
        for easy lookup when altering abundances.
    """

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
        #  Some abstar annotations are blank. Skip these.
        if len(line) < 10:
            len_fail += 1
            continue
        #  Abstar may annotate some sequences as being
        #  kappa or lambda chains. We want heavy chain.
        if 'd_gene' not in line:
            d_fail += 1
            continue
        data_dict = json.loads(line)
        #  If the annotation has no CDR3 sequence,
        #  it's not of much use to us.
        if 'cdr3_nt' not in data_dict:
            cdr3_fail += 1
            continue
        #  Remove any annotations in which the input
        #  sequence has an N in it.
        if data_dict['raw_input'].strip('N').count('N') != 0:
            num_N += 1
            continue
        passed += 1

        #  Save the annotation to the dictionary
        #  using a unique key.
        seq_key = data_dict['seq_id'].split('|')[0]
        if data_dict['prod'] == 'yes':
            annotations['productive'][seq_key] = data_dict
        #  We use out-of-frame sequences exclusively
        #  as our unproductive sequences since these
        #  are least likely to occur from sequencing error
        #  compared to observing a stop codon.
        elif data_dict['junction_in_frame'] == 'no':
            annotations['unproductive'][seq_key] = data_dict
        #  Stop codons in the CDR3 are more likely
        #  occur from vdj recombination and 
        #  junctional diversity than sequencing error alone.
        elif '*' in data_dict['junc_aa']:
            annotations['stop_in_cdr3'][seq_key] = data_dict
        #  Stop codons not in the CDR3 are more likely
        #  to have been from sequencing error. If one
        #  wanted to include them in the error correction
        #  step, it wouldn't be the worst idea.
        else:
            annotations['other'][seq_key] = data_dict

    print('Loading abstar log:'
          '\nLength_fail',len_fail,
          '\nD_gene_fail',d_fail,
          '\nN_filtered', num_N,
          '\nCDR3_fail',cdr3_fail,
          '\nabstar_pass', passed)
    return annotations

def has_shm_indels(annotation: dict) -> bool:
    """Determines if an annotation has shm indels.

    Method to determine whether or not an annotation has SHM indels.
    Any difference between the ungapped and gapped alignments
    is indicative of the presence of SHM indels.

    Parameters
    ----------
    annotation : dict
        Dictionary containing annotation output for a sequence.

    Returns
    -------
    bool
        True if there is any difference between gapped and ungapped alignments.
    """

    return ((annotation['vdj_nt'] != annotation['gapped_vdj_nt'])
       or (annotation['vdj_germ_nt'] != annotation['gapped_vdj_germ_nt']))

def cdr3_position(annotation: dict, junc: bool = False) -> tuple:
    """Determines if reported CDR3 is present in the alignment and returns position if so.

    Abstar is not the most robust of softwares. There are sometimes
    differences between what is in the raw input sequence and what
    abstar reports as the CDR3. Alignment lengths are also not
    completely trustworthy. This method checks whether abstar's
    identified CDR3 was in the alignment of the input and returns the
    start and end positions of the cdr3 in the alignment if so.

    Parameters
    ----------
    annotation : dict
        Dictionary containing annotation output for a sequence.
    junc : bool, optional
        Use junction, otherwise use CDR3 without invariant residues.

    Returns
    -------
    (cdr3_start, cdr3_end) : tuple[str]
        Tuple of start and end positions of CDR3 sequence.
    """

    if junc:
        loc = 'junc_'
    else:
        loc = 'cdr3_'

    try:
        cdr3_start = annotation['vdj_nt'].index(annotation[loc + 'nt'])
        cdr3_end = cdr3_start + len(annotation[loc + 'nt'])
        return cdr3_start,cdr3_end
    except:
        return None

def run_error_correction(annotations: dict, keys: list) -> None:
    """Performs error correction on annotated sequences.

    Error correction merges erroneous sequences by adding their
    abundance to non-erroneous sequences and then removing the errorneous
    sequence. The marginal error correction algorithm seeks to correct
    for sequencing errors that caused large abundance clones to be split
    into many similar sequences. The total error correction algorithm
    targets correction of reverse transcriptase errors. To conserve memory,
    the original dictionary of annotations is altered.

    Parameters
    ----------
    annotation : dict
        Dictionary containing annotation output for a sequence.
    keys : list
        List of string keys in annotations on which error correction
        will be performed. For conservative error correction, perform
        error correction on productives only and then unproductives only.
        For more liberal error correction, you could include sequences
        with stops not in the CDR3, i.e. the key called "other".

    Returns
    -------
    None
    """

    initial_headers = []
    initial_sequences = []
    for key in keys:
        for seq_key in annotations[key]:
            initial_headers.append(annotations[key][seq_key]['seq_id'])
            initial_sequences.append(annotations[key][seq_key]['raw_input'])

    initial_unique_count = len(initial_headers)
    initial_abundance_count = sum([get_abundance(header)
                                   for header in initial_headers])
    marginal_headers = []

    #  Error correction gets performed within a bin of same cprimer, vprimer,
    #  sequence length, time, and replicate.
    grouped_data = group_data(initial_headers, initial_sequences)
    for key in grouped_data:
        marginal_output = error_correct_marginal(grouped_data[key]['headers'],
                                                 grouped_data[key]['sequences'],
                                                 delta_r=1.0,
                                                 delta_a=1.0)
        #  Record headers for logging purposes.
        marginal_headers += marginal_output[0]

        total_output = error_correct_total(marginal_output[0],
                                           marginal_output[1],
                                           d_thresh=2,
                                           a_thresh=None)

        grouped_data[key]['headers'] = total_output[0]
        grouped_data[key]['sequences'] = total_output[1]

    #  Concatenate into single lists for headers and sequences, respectively.
    final_headers = []
    final_sequences = []
    for key in grouped_data:
        final_headers += grouped_data[key]['headers']
        final_sequences += grouped_data[key]['sequences']

    #  Calculate counts at all step of error correction.
    marginal_unique = 0
    marginal_abundance = 0
    for h in marginal_headers:
        abun = get_abundance(h)
        if abun > 0:
            marginal_unique += 1
            marginal_abundance += abun
    final_unique = 0
    final_abundance = 0
    for h in final_headers:
        abun = get_abundance(h)
        if abun > 0:
            final_unique += 1
            final_abundance += abun
    print('-'.join(keys) + '_before_correction', initial_unique_count, initial_abundance_count)
    print('-'.join(keys) + '_after_marginal_correction', marginal_unique, marginal_abundance)
    print('-'.join(keys) + '_after_total_correction', final_unique, final_abundance)

    #  Update abundances and delete 0 abundance sequences.
    for key in keys:
        updated_dict_unique = 0
        updated_dict_abundance = 0
        for h in final_headers:
            seq_key = h.split('|')[0]
            try:
                abun = get_abundance(h)
                if abun > 0:
                    updated_dict_unique += 1
                    updated_dict_abundance += abun
                    annotations[key][seq_key]['seq_id'] = h
                else:
                    del annotations[key][seq_key]
            except:
                pass
        print(key + '_after_type_check', updated_dict_unique, updated_dict_abundance)

def filter_after_error_correction(annotations: dict, annkey: str) -> None:
    """Removes annotations with mismtached primer, shm indels, or corrupt CDR3.

    Parameters
    ----------
    annotations : dict
        Dictionary containing annotation output for a sequence.
    annkey : str
        A key in the annotations dictionary: ['productive', 'unproductive',
        'stop_in_cdr3', 'other'].

    Returns
    -------
    None
    """

    bad_cdr3 = 0
    num_shm_indels = 0
    bad_primer = 0

    for seq_key in list(annotations[annkey].keys()):
        annotation = annotations[annkey][seq_key]
        vprimer = get_vprimer(annotation['seq_id']).split('-')[0]
        v_gene =  annotation['v_gene']['fam']
        if v_gene != vprimer:
            bad_primer += 1
            del annotations[annkey][seq_key]
        elif cdr3_position(annotation) is None:
            bad_cdr3 += 1
            del annotations[annkey][seq_key]
        elif has_shm_indels(annotation):
            num_shm_indels += 1
            del annotations[annkey][seq_key]

    print(annkey + '_bad_cdr3', bad_cdr3)
    print(annkey + '_num_shm_indels', num_shm_indels)
    print(annkey + '_bad_primer', bad_primer)

def annotations_pipeline(infile: str) -> dict:
    """Loads abstar file and saves annotations after error correction and filtering.

    Parameters
    ----------
    infile : (str)
        Path to abstar file.

    Returns
    -------
    None
    """

    patient = infile.split('/')[-1].split('.')[0]
    annotations = load_abstar(infile)
    for key in annotations:
        abun = sum([get_abundance(annotations[key][seq_key]['seq_id'])
                    for seq_key in annotations[key]])
        print('initial_' + key, len(annotations[key]), abun)

    run_error_correction(annotations, ['productive'])
    run_error_correction(annotations, ['unproductive'])

    del annotations['stop_in_cdr3']
    del annotations['other']

    filter_after_error_correction(annotations, 'productive')
    filter_after_error_correction(annotations, 'unproductive')
    for key in annotations:
        abun = sum([get_abundance(annotations[key][seq_key]['seq_id'])
                    for seq_key in annotations[key]])
        print('final_' + key, len(annotations[key]), abun)

    for key in annotations:
        annotations[key] = list(annotations[key].values())

    return annotations

def make_lineages(infile: str) -> dict:
    """Loads annotations file and creates lineages by performing SLC on VJL bins.

    Parameters
    ----------
    infile : str
        Path to annotations file.

    Returns
    -------
    lineages : dict
        Nested dictionaries of clustered annotations [V][J][L][cluster_id].
    """

    patient = infile.split('/')[-1].split('_')[0]
    annotations = json_open(infile)
    productive_lineages = vjl_slc(annotations['productive'],abstar=True)
    unproductive_lineages = vjl_slc(annotations['unproductive'],abstar=True)
    lineages = {'productive': productive_lineages,
                'unproductive': unproductive_lineages}
    return lineages

def update_header(header: str, abundance: int) -> str:
    """Changes abundance count in header.

    Parameters
    ----------
    header : str
        String containing information about a sequence/annotation,
    abundance : int
        Abundance count to replace that currently in header.

    Returns
    -------
    updated_header : str
        String containing updated abundance.
    """

    header_split = header.split('|')
    header_split[3] = 'DUPCOUNT=' + str(abundance)
    updated_header = '|'.join(header_split)
    return updated_header

def get_naive_cdr3(annotation: dict,  nt: bool = True, junc: bool = True) -> str:
    """Obtains naive CDR3 from abstar annotation.

    Parameters
    ----------
    annotation : dict
        Dictionary containing annotation output for a sequence.
    nt : bool, optional
        Use nucleotide output, otherwise amino acid.
    junc : bool, optional
        Use junction, otherwise use CDR3 without invariant residues.

    Returns
    -------
    naive_cdr3 : str
        CDR3 of naive sequence.
    """

    if junc:
        loc = 'junc_'
    else:
        loc = 'cdr3_'

    cdr3_positions = cdr3_position(annotation, junc=junc)
    if cdr3_positions is not None:
        if nt == True:
            return annotation['vdj_germ_nt'][cdr3_positions[0]:cdr3_positions[1]]
        else:
            return translate(annotation['vdj_germ_nt'][cdr3_positions[0]:cdr3_positions[1]])
    else:
       return None

def merge_replicates_within_lineage(lineage: list) -> list:
    """Merges identical sequences from different replicates at the same timepoint for a single lineage.

    Parameters
    ----------
    lineage : list
        List of annotations.

    Returns
    -------
    condensed_lineage : list
        List of annotations after merging replicates.
    """

    #  Create dictionary to identify identical sequences across replicates
    #  at the same timepoint.
    replicate_dict = {}
    for annotation in lineage:
        sequence = annotation['raw_input']
        time = get_time(annotation['seq_id'])
        replicate_key = (sequence, time)
        if replicate_key not in replicate_dict:
            replicate_dict[replicate_key] = []
        replicate_dict[replicate_key].append(annotation)

    condensed_lineage = []
    for key in replicate_dict:
        if len(replicate_dict[key]) > 1:
            abundance = sum([get_abundance(annotation['seq_id'])
                             for annotation in replicate_dict[key]])
            annotation_0 = replicate_dict[key][0]
            naive_cdr3_0 = get_naive_cdr3(annotation_0)
            replicate_0 = get_replicate(annotation_0['seq_id'])

            #  Check that the sequence produces the same annotation results
            #  across replicates.
            for annotation in replicate_dict[key][1:]:
                naive_cdr3 = get_naive_cdr3(annotation)
                if naive_cdr3 != naive_cdr3_0:
                    replicate = get_replicate(annotation['seq_id'])
                    print('Different annotation results for these replicates:', replicate_0, replicate)
                    print(replicate_0, annotation_0['raw_input'], naive_cdr3_0)
                    print(replicate, annotation['raw_input'], naive_cdr3)
                    print('\n\n')
            annotation_0['seq_id'] = update_header(annotation_0['seq_id'], abundance)
            condensed_lineage.append(annotation_0)
        else:
            condensed_lineage.append(replicate_dict[key][0])
    return condensed_lineage

def merge_replicates(lineages: dict, productive: bool = True) -> dict:
    """Apply replicate merging to each lineage.

    Parameters
    ----------
    lineage : dict
        Nested dictionaries of lineages [V][J][L][cluster_id].
    productive : bool, optional
        Bool to specify whether to keep lineages with unproductive progenitors.
        Note: productive lineages can have unproductive progenitors due to uncertainty.

    Returns
    -------
    condensed_lineages : dict
        Dictionary of lineages by [(V, J, L, cluster_id)].
    """

    condensed_lineages = {}
    if productive:
        for v in lineages:
            for j in lineages[v]:
                for l in lineages[v][j]:
                    for cluster_id in lineages[v][j][l]:
                        merged_lineage = merge_replicates_within_lineage(lineages[v][j][l][cluster_id])
                        if get_lineage_progenitor_cdr3(merged_lineage) == '':
                            continue
                        condensed_lineages[(v,j,l,cluster_id)] = merge_replicates_within_lineage(lineages[v][j][l][cluster_id])
    else:
        for v in lineages:
            for j in lineages[v]:
                for l in lineages[v][j]:
                    for cluster_id in lineages[v][j][l]:
                        condensed_lineages[(v,j,l,cluster_id)] = merge_replicates_within_lineage(lineages[v][j][l][cluster_id])
    return condensed_lineages

def denest_lineages(lineages: dict, productive: bool = True) -> dict:
    """Convert lineage dictionary to one that doesn't nest for simpler looping.

    Parameters
    ----------
    lineage : dict
        Nested dictionaries of lineages [V][J][L][cluster_id].
    productive : bool, optional
        Bool to specify whether to keep lineages with unproductive progenitors.
        Note: productive lineages can have unproductive progenitors due to uncertainty.

    Returns
    -------
    condensed_lineages : dict
        Dictionary of lineages by [(V, J, L, cluster_id)].
    """

    condensed_lineages = {}
    if productive:
        for v in lineages:
            for j in lineages[v]:
                for l in lineages[v][j]:
                    for cluster_id in lineages[v][j][l]:
                        if get_lineage_progenitor_cdr3(lineages[v][j][l][cluster_id]) == '':
                            continue
                        condensed_lineages[(v,j,l,cluster_id)] = lineages[v][j][l][cluster_id]
    else:
        for v in lineages:
            for j in lineages[v]:
                for l in lineages[v][j]:
                    for cluster_id in lineages[v][j][l]:
                        condensed_lineages[(v,j,l,cluster_id)] = lineages[v][j][l][cluster_id]
    return condensed_lineages

def get_lineage_progenitor(lineage: list) -> str:
    """Obtains the most common VDJ germline sequence.

    Parameters
    ----------
    lineage : list
        List of annotations.

    Returns
    -------
    most_common_naive : str
        Most common VDJ germline sequence in a lineage.

    """

    germline_counter = Counter([annotation['vdj_germ_nt'] for annotation in lineage])
    most_common_naive = germline_counter.most_common()[0][0]
    return most_common_naive

def get_lineage_progenitor_cdr3(lineage: list, nt: bool = True, junc: bool = True) -> str:
    """Obtains the most common naive cdr3 not containg a stop codon from a lineage.

    If the most common naive sequence of a lineage contains a stop
    codon, the progenitor of the lineage is chosen iteratively
    by examining the next most common naive CDR3 until it does not
    contain any stop codons. Otherwise, an empty string is returned
    indicating all naive CDR3s in this lineage contain a stop codon.

    Parameters
    ----------
    lineage : list
        List of annotations.
    nt : bool, optional
        Use nucleotide output, otherwise amino acid.
    junc : bool, optional
        Use junction, otherwise use CDR3 without invariant residues.

    Returns
    -------
    str
        Most common naive CDR3 in a lineage. If no such productive CDR3 exists,
        returns an empty string instead.
    """

    naive_cdr3s = []
    for annotation in lineage:
        naive_cdr3s.append(get_naive_cdr3(annotation, nt=nt, junc=junc))
    counter = sort_dict_by_value(Counter(naive_cdr3s))
    for key in counter:
        if key is None:
            return ''
        if nt:
            if '*' in translate(key):
                continue
            else:
                return key
        else:
            if '*' in key:
                continue
            else:
                return key
    return ''

def create_sonia_input(infile: str) -> pd.DataFrame:
    """Generates a csv of information about lineage progenitors to use as input for SONIA.

    Lineages first undergo replicate merging. This is done so that only completely
    unique sequences are used at each timepoint when determining the most common
    naive CDR3. Once the most common CDR3 is obtained, it's inspected
    to ensure it either doesn't contain any stop codons and that is begins
    with a cysteine (C) and ends in a tryptophan (W) or phenylalanine (F).
    SONIA works on naive productive sequences only and the invariant residues
    are necessary due to SONIA's preinstalled CDR3 anchors. Finally, only
    unique nucleotide rearrangements are kept. (I.e. multiple combinations
    of (CDR3, V, and J) may be in the file, but know that they ensue from
    different nucleotide rearrangement events.)

    Parameters
    ----------
    infile : str
        Path to file containing lineages.

    Returns
    -------
    sonia_df : pandas.DataFrame
        DataFrame object containing information for SONIA.
    """

    in_lineages = json_open(infile)['productive']
    lineages = merge_replicates(in_lineages)
    sonia_input = []
    patient = infile.split('/')[-1].split('_')[0]

    for index,key in enumerate(lineages):
        progenitor_cdr3 = get_lineage_progenitor_cdr3(lineages[key])
        if progenitor_cdr3 is None:
            continue
        if len(progenitor_cdr3) == 0:
            continue
        progenitor_cdr3_aa = translate(progenitor_cdr3)
        if ((progenitor_cdr3_aa[0] == 'C')
           and (progenitor_cdr3_aa[-1] == 'W' or progenitor_cdr3_aa[-1] == 'F')):
            sonia_input.append((progenitor_cdr3_aa, key[0], key[1],
                                len(lineages[key]), patient, progenitor_cdr3))
    df_sonia = pd.DataFrame(sonia_input, columns=['progenitor_cdr3', 'v_gene',
                                              'j_gene', 'lineage_size',
                                              'patient', 'nt_cdr3'])
    #  Sort by lineage size so smaller duplicate nucleotide recombination
    #  events get removed.
    df_sonia = df_sonia.sort_values(by=['lineage_size'],ascending=False)

    #  Keep independent nucleotide recombination events only
    df_sonia.drop_duplicates(subset=['v_gene', 'j_gene', 'progenitor_cdr3', 'nt_cdr3'],
                             keep='first')
    return df_sonia

def main():
    import argparse

    parser = argparse.ArgumentParser(
        description='Pipeline to process abstar output:'
                    'error correction on sequences and filter annotations;'
                    'make lineages, merge replicates in lineages,'
                    'get progenitor CDR3s; create csv to be used as input'
                    'for SONIA and expansion analysis.')
    parser.add_argument('--annotations', dest='annotations_input', nargs=2,
                        metavar=('PATH/TO/FOLDER/FILE', 'PATH/TO/OUTFILE'),
                        help='path to abstar output (.json), '
                        'path for annotations output (.json)')
    parser.add_argument('--lineages', dest='lineages_input', nargs=2,
                        metavar=('PATH/TO/FOLDER/FILE', 'PATH/TO/OUTFILE'),
                        help='path to annotations output (.json), '
                        'path for lineages output (.json)')
    parser.add_argument('--sonia', dest='sonia_input', nargs=2,
                        metavar=('PATH/TO/FOLDER/FILE', 'PATH/TO/OUTFILE'),
                        help='path to lineages file (.json), '
                        'path for SONIA input (.csv)')
    args = parser.parse_args()

    if args.annotations_input is not None:
        annotations = annotations_pipeline(args.annotations_input[0])
        json_save(args.annotations_input[1], annotations)
    if args.lineages_input is not None:
        lineages = make_lineages(args.lineages_input[0])
        json_save(args.lineages_input[1], lineages)
    if args.sonia_input is not None:
        df_sonia = create_sonia_input(args.sonia_input[0])
        df_sonia.to_csv(args.sonia_input[1], index=False)

if __name__=='__main__':
    main()

