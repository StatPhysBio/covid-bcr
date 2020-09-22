#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Script to combine all samples for a patient into one FASTA file with annotated times, severitie, replicates, and patient ID.
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

from utils import *

def annotate_fasta(fasta: str, singletons: bool = True) -> (list, list):
    """Annotates a online deduplicated FASTA file from the pRESTO output with patient, replicate, and time information.

    Replicate number will be taken from the filename. If no replicate number is supplied in the filename
    (the number after the dash), the replicate number will be set to 0.

    Parameters
    ----------
    fasta : str
        Path to FASTA.
    singletons : bool, optional
        Specifies whether or not to keep singletons.

    Returns
    -------
    sequences : list
        List of sequences from FASTA file.
    headers : list
        List of annotated headers from FASTA file.
    """

    sequences = []
    headers = []

    filename = fasta.split("/")[-1]
    sample = 'IgG' + filename.split("_")[0].replace("S","")

    for seq in parse(fasta, 'fasta'):
        uid, cprimer, vprimer, abundance = (x for x in seq.id.split('|'))

        abundance_int = int(abundance.split("=")[-1])
        if not singletons:
            if abundance_int == 1:
                continue

        header_info = [uid+"="+CONST_SAMPLE_DICT[sample]['patient'],
                       cprimer, vprimer,abundance,
                       "TIME=" + str(CONST_SAMPLE_DICT[sample]['sample day'][0]),
                       "SEVERITY=" + str(CONST_SAMPLE_DICT[sample]['severity']),
                       "REPLICATE=" + str(CONST_SAMPLE_DICT[sample]['replicate'][0])]

        header = "|".join(header_info)
        header = ">" + header
        headers.append(header)
        sequences.append(str(seq.seq))
    return sequences, headers

def map_sample_to_patient(sample_files: list) -> dict:
    """Groups all sample files by patient.

    Parameters
    ----------
    sample_files : list
        List of pRESTO pre-processed files.

    Returns
    -------
    files_by_patient : dict
        Dictionary with keys given by patient ID and values which are sample files.
    """
    files_by_patient = {}
    for patient in CONST_DATA_DICT:
        files_by_patient[patient] = {}
        for f in sample_files:
            sample = f.split('/')[-1].split('_')[0].replace('S', '')
            sample_w_igg = 'IgG' + sample
            print(sample)
            if sample_w_igg in CONST_DATA_DICT[patient]['sample']:
                files_by_patient[patient][sample] = f
    files_by_patient = sort_dict(files_by_patient)
    for key in files_by_patient:
        files_by_patient[key] = sort_dict(files_by_patient[key])
    return files_by_patient

def assemble_combined_fasta(savename: str, sample_files: list) -> None:
    """Annotates given files and saves them all to a single FASTA file.

    Parameters
    ----------
    savename : str
        Path and name of FASTA file to save.
    sample_files : list
        List of samples files to annotate and write to a single FASTA file.

    Returns
    -------
    None
    """

    all_sequences = []
    all_headers = []
    for f in sample_files:
        sequences, headers = annotate_fasta(f)
        all_sequences += sequences
        all_headers += headers
    create_new_fasta(savename, all_headers, all_sequences)

def main():
    import argparse

    parser = argparse.ArgumentParser(
        description='Combine multiple patient samples into a single FASTA file.')
    parser.add_argument('--dirs', nargs='+', type=str,
                        help='directories which contain pRRESTO collapsed files')
    parser.add_argument('--savedir', type=str, help='directory to save files')
    parser.add_argument('--bcellinfo', type=str, default='total_b_cell_info.csv',
                        help='path to SRA-style csv file')
    args = parser.parse_args()

    dirs = args.dirs
    savedir = args.savedir
    global CONST_DATA_DICT
    CONST_DATA_DICT = get_bcell_info(args.bcellinfo)
    global CONST_SAMPLE_DICT
    CONST_SAMPLE_DICT = get_bcell_info(args.bcellinfo, group='sample_name')
    if dirs is None:
        print("No directory specified, so no files will be located.")
        return
    if savedir is None:
        print("No save directory is specified.")
        return


    files = []
    for d in dirs:
        files += get_files(d,"unique.fasta")
    files_by_patient  = map_sample_to_patient(files)
    print(files_by_patient)
    for patient in files_by_patient:
        print(files_by_patient[patient])
        assemble_combined_fasta(savedir + patient + '.fasta', list(files_by_patient[patient].values()))

if __name__ == '__main__':
    main()
