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

#  Dictionary used to map samples to patients when assembling FASTAs together.
CONST_SAMPLE_DICT = {'1': {'IgG-nCoV-RBD': -999.0, 'IgM-nCoV-RBD': -999.0,
                           'sample day': 0, 'round': 0, 'severity': 'Healthy', 'patient': 'H1'},
                     '2': {'IgG-nCoV-RBD': -999.0, 'IgM-nCoV-RBD': -999.0,
                           'sample day': 0, 'round': 0, 'severity': 'Healthy', 'patient': 'H2'},
                     '3': {'IgG-nCoV-RBD': -999.0, 'IgM-nCoV-RBD': -999.0,
                           'sample day': 0, 'round': 0, 'severity': 'Healthy', 'patient': 'H3'},
                     '4': {'IgG-nCoV-RBD': 0.33, 'IgM-nCoV-RBD': 0.1315,
                           'sample day': 6, 'round': 1, 'severity': 'Severe', 'patient': '19'},
                     '5': {'IgG-nCoV-RBD': 0.2745, 'IgM-nCoV-RBD': 0.361,
                           'sample day': 8, 'round': 1, 'severity': 'Moderate', 'patient': '3'},
                     '6': {'IgG-nCoV-RBD': 0.799, 'IgM-nCoV-RBD': 0.8005,
                           'sample day': 2, 'round': 1, 'severity': 'Mild', 'patient': '2'},
                     '7': {'IgG-nCoV-RBD': 1.477, 'IgM-nCoV-RBD': 1.591,
                           'sample day': 10, 'round': 1, 'severity': 'Moderate', 'patient': '5'},
                     '8': {'IgG-nCoV-RBD': 1.89, 'IgM-nCoV-RBD': 2.904,
                           'sample day': 13, 'round': 1, 'severity': 'Moderate', 'patient': '6'},
                     '9': {'IgG-nCoV-RBD': 1.447, 'IgM-nCoV-RBD': 1.9265,
                           'sample day': 15, 'round': 1, 'severity': 'Mild', 'patient': '2'},
                     '10': {'IgG-nCoV-RBD': 0.1645, 'IgM-nCoV-RBD': 0.3355,
                            'sample day': 11, 'round': 1, 'severity': 'Moderate', 'patient': '7'},
                     '11': {'IgG-nCoV-RBD': 1.218, 'IgM-nCoV-RBD': 2.3245,
                            'sample day': 16, 'round': 1, 'severity': 'Moderate', 'patient': '7'},
                     '12': {'IgG-nCoV-RBD': 2.0025, 'IgM-nCoV-RBD': 1.139,
                            'sample day': 8, 'round': 1, 'severity': 'Severe', 'patient': '18'},
                     '13': {'IgG-nCoV-RBD': 2.028, 'IgM-nCoV-RBD': 2.1485,
                            'sample day': 27, 'round': 1, 'severity': 'Moderate', 'patient': '5'},
                     '14': {'IgG-nCoV-RBD': 2.08, 'IgM-nCoV-RBD': 1.7305,
                            'sample day': 38, 'round': 1, 'severity': 'Moderate', 'patient': '3'},
                     '15': {'IgG-nCoV-RBD': 1.4265, 'IgM-nCoV-RBD': 2.0875,
                            'sample day': 34, 'round': 1, 'severity': 'Mild', 'patient': '2'},
                     '16': {'IgG-nCoV-RBD': 2.293, 'IgM-nCoV-RBD': 2.7875,
                            'sample day': 28, 'round': 1, 'severity': 'Moderate', 'patient': '6'},
                     '17': {'IgG-nCoV-RBD': 1.8605, 'IgM-nCoV-RBD': 2.2155,
                            'sample day': 30, 'round': 1, 'severity': 'Severe', 'patient': '18'},
                     '18': {'IgG-nCoV-RBD': 2.092, 'IgM-nCoV-RBD': 3.0465,
                            'sample day': 39, 'round': 1, 'severity': 'Moderate', 'patient': '7'},
                     '19': {'IgG-nCoV-RBD': 1.892, 'IgM-nCoV-RBD': 1.8885,
                            'sample day': 18, 'round': 2, 'severity': 'Moderate', 'patient': '4'},
                     '20': {'IgG-nCoV-RBD': 2.3095, 'IgM-nCoV-RBD': 2.129,
                            'sample day': 14, 'round': 2, 'severity': 'Moderate', 'patient': '8'},
                     '21': {'IgG-nCoV-RBD': 1.982, 'IgM-nCoV-RBD': 1.823,
                            'sample day': 22, 'round': 2, 'severity': 'Mild', 'patient': '1'},
                     '22': {'IgG-nCoV-RBD': 1.9595, 'IgM-nCoV-RBD': 2.7705,
                            'sample day': 32, 'round': 2, 'severity': 'Moderate', 'patient': '8'},
                     '23': {'IgG-nCoV-RBD': 2.276, 'IgM-nCoV-RBD': 2.4285,
                            'sample day': 34, 'round': 2, 'severity': 'Moderate', 'patient': '4'},
                     '24': {'IgG-nCoV-RBD': 2.0465, 'IgM-nCoV-RBD': 1.581,
                            'sample day': 43, 'round': 2, 'severity': 'Mild', 'patient': '1'},
                     '25': {'IgG-nCoV-RBD': 0.4795, 'IgM-nCoV-RBD': 0.904,
                            'sample day': 5, 'round': 2, 'severity': 'Moderate', 'patient': '9'},
                     '26': {'IgG-nCoV-RBD': 0.1525, 'IgM-nCoV-RBD': 0.2125,
                            'sample day': 8, 'round': 2, 'severity': 'Moderate', 'patient': '9'},
                     '27': {'IgG-nCoV-RBD': 2.023, 'IgM-nCoV-RBD': 2.196,
                            'sample day': 37, 'round': 2, 'severity': 'Moderate', 'patient': '8'},
                     '28': {'IgG-nCoV-RBD': 1.551, 'IgM-nCoV-RBD': 2.074,
                            'sample day': 18, 'round': 2, 'severity': 'Moderate', 'patient': '9'},
                     '29': {'IgG-nCoV-RBD': 1.97, 'IgM-nCoV-RBD': 2.107,
                            'sample day': 44, 'round': 3, 'severity': 'Severe', 'patient': '16'},
                     '30': {'IgG-nCoV-RBD': 2.034, 'IgM-nCoV-RBD': 1.254,
                            'sample day': 53, 'round': 3, 'severity': 'Severe', 'patient': '16'},
                     '31': {'IgG-nCoV-RBD': 2.4255, 'IgM-nCoV-RBD': 2.848,
                            'sample day': 22, 'round': 3, 'severity': 'Severe', 'patient': '17'},
                     '32': {'IgG-nCoV-RBD': 2.029, 'IgM-nCoV-RBD': 2.886,
                            'sample day': 36, 'round': 3, 'severity': 'Severe', 'patient': '17'},
                     '33': {'IgG-nCoV-RBD': 0.159, 'IgM-nCoV-RBD': 1.311,
                            'sample day': 7, 'round': 3, 'severity': 'Moderate', 'patient': '12'},
                     '34': {'IgG-nCoV-RBD': 0.152, 'IgM-nCoV-RBD': 1.1875,
                            'sample day': 11, 'round': 3, 'severity': 'Moderate', 'patient': '12'},
                     '35': {'IgG-nCoV-RBD': 0.123, 'IgM-nCoV-RBD': 0.1975,
                            'sample day': 7, 'round': 3, 'severity': 'Moderate', 'patient': '13'},
                     '36': {'IgG-nCoV-RBD': 1.594, 'IgM-nCoV-RBD': 1.0525,
                            'sample day': 13, 'round': 3, 'severity': 'Moderate', 'patient': '13'},
                     '37': {'IgG-nCoV-RBD': 0.185, 'IgM-nCoV-RBD': 0.178,
                            'sample day': 3, 'round': 3, 'severity': 'Moderate', 'patient': '14'},
                     '38': {'IgG-nCoV-RBD': 0.308, 'IgM-nCoV-RBD': 0.687,
                            'sample day': 7, 'round': 3, 'severity': 'Moderate', 'patient': '14'},
                     '39': {'IgG-nCoV-RBD': 0.248, 'IgM-nCoV-RBD': 0.3735,
                            'sample day': 8, 'round': 3, 'severity': 'Moderate', 'patient': '10'},
                     '40': {'IgG-nCoV-RBD': 0.7375, 'IgM-nCoV-RBD': 0.726,
                            'sample day': 14, 'round': 3, 'severity': 'Moderate', 'patient': '10'},
                     '41': {'IgG-nCoV-RBD': 2.0945, 'IgM-nCoV-RBD': 1.0145,
                            'sample day': 39, 'round': 3, 'severity': 'Severe', 'patient': '15'},
                     '42': {'IgG-nCoV-RBD': 1.9905, 'IgM-nCoV-RBD': 0.8055,
                            'sample day': 41, 'round': 3, 'severity': 'Severe', 'patient': '15'},
                     '43': {'IgG-nCoV-RBD': 2.0735, 'IgM-nCoV-RBD': 0.6615,
                            'sample day': 71, 'round': 3, 'severity': 'Severe', 'patient': '15'},
                     '44': {'IgG-nCoV-RBD': 0.1635, 'IgM-nCoV-RBD': 0.4795,
                            'sample day': 10, 'round': 3, 'severity': 'Moderate', 'patient': '11'},
                     '45': {'IgG-nCoV-RBD': 0.518, 'IgM-nCoV-RBD': 1.5075,
                            'sample day': 15, 'round': 3, 'severity': 'Moderate', 'patient': '11'}}

def annotate_fasta(fasta: str, patient: str = '', timepoint: str = '',
                   severity: str = '', replicate: str = '',
                   oneline_collapsed: bool = False,
                   singletons: bool = True) -> (list, list):
    """Annotates a online deduplicated FASTA file from the pRESTO output with patient, replicate, and time information.

    Replicate number will be taken from the filename. If no replicate number is supplied in the filename
    (the number after the dash), the replicate number will be set to 0.

    Parameters
    ----------
    fasta : str
        Path to FASTA.
    patient : str, optional
        Patient ID.
    timepoint : str, optional
        Days after infection of an individual.
    severity : str, optional
        Severity of the disease in the individual.
    replicate : str, optional
        Replicate number.
    oneline_collapsed : bool, optional
        Specifies if this is output from pRESTO.
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

    if oneline_collapsed and replicate == '':
        filename = fasta.split("/")[-1]
        sample = filename.split("_")[0].replace("S","")
        if "-" in sample:
            replicate = sample.split("-")[-1]
        else:
            replicate = '0'

    for seq in parse(fasta, 'fasta'):
        uid, cprimer, vprimer, abundance = (x for x in seq.id.split('|'))

        abundance_int = int(abundance.split("=")[-1])
        if not singletons:
            if abundance_int == 1:
                continue

        if oneline_collapsed and not replicate:
            sample = sample.split("-")[0]
            header_info = [uid+"="+CONST_SAMPLE_DICT[sample]['patient'],
                           cprimer, vprimer,abundance,
                           "TIME=" + str(CONST_SAMPLE_DICT[sample]['sample day']),
                           "SEVERITY=" + str(CONST_SAMPLE_DICT[sample]['severity']),
                           "REPLICATE=" + replicate]
        if oneline_collapsed and replicate:
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
            sample_and_replicate = f.split("/")[-1].split("_")[0]
            sample = sample_and_replicate.split("-")[0].replace("S","")
            if sample in CONST_DATA_DICT[patient]['sample']:
                files_by_patient[patient][sample_and_replicate] = f
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
        sequences, headers = annotate_fasta(f, oneline_collapsed=True)
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
    args = parser.parse_args()

    dirs = args.dirs
    savedir = args.savedir

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

    for patient in files_by_patient:
        assemble_combined_fasta(savedir + patient + '.fasta', list(files_by_patient[patient].values()))

if __name__ == '__main__':
    main()
