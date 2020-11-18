#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Script to assemble FASTA file for input to IGoR.
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

from abstar_pipeline import (merge_replicates, get_lineage_progenitor,
                             get_lineage_progenitor_cdr3, translate)
from utils import *

def assemble_fasta_for_igor(infiles: list, outfile: str) -> None:
    """Assembles all unproductive lineage progenitors into one FASTA for IGoR input.

    CDR3s are constrained to begin with a C and end in a W.

    Parameters
    ----------
    infiles : list
        List of paths to all lineages.
    outfile : str
        Path of FASTA to be saved.

    Returns
    -------
    None
    """

    sequences_for_igor = []
    headers_for_igor = []
    for f in infiles:
        lineages = json_open(f)['unproductive']
        lineages = merge_replicates(lineages, productive=False)
        patient = get_patient(lineages[list(lineages.keys())[0]][0]['seq_id'])
        for lin in lineages:
            lineage_progenitor = get_lineage_progenitor(lineages[lin])
            lineage_progenitor_cdr3 = get_lineage_progenitor_cdr3(lineages[lin])
            first_aa = translate(lineage_progenitor_cdr3[0:3])
            last_aa = translate(lineage_progenitor_cdr3[-3:])
            if first_aa != 'C' or last_aa != 'W':
                continue
            sequences_for_igor.append(lineage_progenitor)
            header = "-".join(list(lin)+[patient])
            headers_for_igor.append(header)

    create_new_fasta(outfile, headers_for_igor, sequences_for_igor)

def main():
    import argparse

    parser = argparse.ArgumentParser(
        description='Combine all unproductive lineage progenitors to make IGOR input.')
    parser.add_argument('--infiles', nargs='+', type=str,
                        help='.json lineage files. e.g. usage: /PATH/TO/LINEAGES/*')
    parser.add_argument('--outfile', type=str,
                       help='path and name of FASTA file for IGoR input.')
    args = parser.parse_args()

    infiles = args.infiles
    outfile = args.outfile

    assemble_fasta_for_igor(infiles, outfile)

if __name__ == '__main__':
    main()
