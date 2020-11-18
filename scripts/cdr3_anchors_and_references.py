#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Script to make CDR3 anchor files and genomic files for IGoR and SONIA.
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

from collections import OrderedDict
from operator import itemgetter

from Bio.SeqIO import parse

def sort_dict(unsorted_dict: dict) -> OrderedDict:
    """Sorts a dictionary by keys.

    Parameters
    ----------
    unsorted_dict : dict

    Returns
    -------
    OrderedDict
        Dictionary sorted by keys.
    """

    return OrderedDict(sorted(unsorted_dict.items(), key=itemgetter(0)))

def get_v_gaps(header: str) -> int:
    """Returns number of gaps in header of IMGT gapped V gene FASTA.

    Parameters
    ----------
    header : str
        String of information for a single V gene from IMGT gapped V gene FASTA.

    Returns
    -------
    int
        Gaps in header.
    """

    col = header[12]
    return int(col[col.find('+') + 1:col.find('=')])

def get_v_anchor(gaps: int) -> str:
    """Returns V gene CDR3 anchor position.

    Parameters
    ----------
    gaps : int
        Number of gaps from a header in IMGT gapped V gene FASTA.

    Returns
    -------
    str
        V gene CDR3 anchor position.
    """

    return str(3 * (104 - 1) - gaps)

def get_j_gaps(header: str) -> int:
    """Returns number of gaps in header of IMGT gapped J gene FASTA.

    Parameters
    ----------
    header : str
        String of information for a single V gene from IMGT gapped J gene FASTA.

    Returns
    -------
    int
        Gaps in header.
    """

    col = header[12]
    gaps = int(col.split('+')[0])
    return gaps

def get_j_anchor(gaps: int) -> str:
    """Returns J gene CDR3 anchor position.

    Parameters
    ----------
    gaps : int
        Number of gaps from a header in IMGT gapped V gene FASTA.

    Returns
    -------
    str
        V gene CDR3 anchor position.
    """

    return str(gaps - 34)

def get_function(header: str) -> str:
    """Returns the type of functionality of a gene.

    Parameters
    ----------
    header : str
        String of information for a single V gene from IMGT gapped gene FASTA.

    Returns
    -------
    str
        Type of functionality.
    """

    return header[3]

def trim_dict(gene_dict: dict) -> dict:
    """Saves only the first two alleles in a gene.

    There's no need for all the alleles from a gene in order for IGoR
    to work. In fact, using all the alleles might make IGoR take longer
    without any improvement in identification since we're not so
    concerned with the specific allele as much as the gene.

    Parameters
    ----------
    gene_dict : dict
        Nested dictionary of gene and then allele keys which contains information
        about the gaps and functionality for each allele.

    Returns
    -------
    out_dict : dict
        Dictionary which has a single key for gene+allele combinations.
        It contains  only two alleles for each gene.
    """

    gene_dict = sort_dict(gene_dict)
    out_dict = {}
    for gene in gene_dict:
        gene_dict[gene] = sort_dict(gene_dict[gene])
        for idx, allele in enumerate(gene_dict[gene]):
            if idx > 1:
                break
            out_dict[gene + '*' + allele] = gene_dict[gene][allele]
    return out_dict

def get_gene_dict(gapped_fasta: str, gene: str = 'V') -> dict:
    """Creates a dictionary of information taken from a gapped FASTA file.

    Parameters
    ----------
    gapped_fasta : str
        Path to FASTA file with gaps.
    gene : str, optional
        Specifies which gene to use. Options = ['V', 'J']. Default = 'V'.

    Returns
    -------
    gene_dict : dict
        Dictionary which has a single key for gene+allele combinations.
        It contains only two alleles for each gene. For each gene+allele,
        it contains the numbers of gaps and type of functionality.
    """

    gene_dict = {}
    if gene == 'V':
        gap_func = get_v_gaps
        anchor_func = get_v_anchor
    elif gene == 'J':
        gap_func = get_j_gaps
        anchor_func = get_j_anchor

    with open(gapped_fasta) as infile:
        for line in infile:
            if ">" in line:
                line = line.split('|')
                gene = line[1].split('*')[0]
                if gene not in gene_dict:
                    gene_dict[gene] = {}
                allele = line[1].split('*')[-1]
                gaps = gap_func(line)
                gene_dict[gene][allele] = {'gaps':anchor_func(gaps),
                                           'func':get_function(line)}
    return trim_dict(gene_dict)

def write_to_csv(outfile: str, gene_dict: dict, program: str = 'igor') -> None:
    """Writes gene CDR3 anchor information to a csv for a given software input.

    IGoR uses ; as a delimiter. SONIA uses , as a delimiter.

    Parameters
    ----------
    outfile : str
        Path for the csv file to be written.
    gene_dict : dict
        Dictionary which has a single key for gene+allele combinations.
        It contains only two alleles for each gene. For each gene+allele,
        it contains the numbers of gaps and type of functionality.
    program : str, optional
        Specifies either IGoR or SONIA. This specifies delimiter and output information.
        Options = ['igor', 'sonia']. Default = 'igor'

    Returns
    -------
    None
    """

    if program.lower() == 'igor':
        delimiter = ';'
        with open(outfile, 'w') as outfile:
            outfile.write('gene' + delimiter + 'anchor_index\n')
            for key in gene_dict:
                outfile.write(key + delimiter + gene_dict[key]['gaps'] + "\n")
    elif program.lower() == 'sonia':
        delimiter = ','
        with open(outfile, 'w') as outfile:
            outfile.write('gene' + delimiter
                          + 'anchor_index' + delimiter
                          + 'function\n')
            for key in gene_dict:
                outfile.write(key + delimiter
                              + gene_dict[key]['gaps'] + delimiter
                              + gene_dict[key]['func'] + "\n")
    else:
        print("Option not recognized. No file written.")

def get_genomic_ref(fasta: str, outfile: str, gene_dict:dict) -> None:
    """Parses ungapped FASTA file and writes genes+alleles to FASTA file.

    This function uses only the alleles which were selected from the gapped FASTA.
    This ensures that the anchor file and genomic reference file correspond
    to one another.

    Parameters
    ----------
    fasta : str
        Ungapped FASTA file.
    outfile : str
       Path for the FASTA file to be written.
    gene_dict : dict
        Dictionary which has a single key for gene+allele combinations.
        It contains only two alleles for each gene. For each gene+allele,
        it contains the numbers of gaps and type of functionality.

    Returns
    -------
    None
    """

    seqs = []
    headers = []
    for seq in parse(fasta, 'fasta'):
        headers.append(seq.description)
        seqs.append(str(seq.seq))

    with open(outfile, 'w') as f:
        for idx,header in enumerate(headers):
            if header in gene_dict:
                f.write(">"+header+"\n")
                f.write(seqs[idx] + "\n")

def main():
    import argparse

    parser = argparse.ArgumentParser(
        description='creates genomic references and CDR3 anchors '
        'from abstar references. Files must be named specifically:'
        'abstar_gapped_V.fasta, abstar_gapped_J.fasta,'
        'abstar_genomic_Vs.fasta, abstar_genomic_Ds.fasta, abstar_genomic_Js.fasta')
    parser.add_argument('--indir', type=str, help='path to reference files.'
                        ' This will also be where the output is saved.')
    args = parser.parse_args()

    genes = ['V', 'J']
    softwares = ['igor', 'sonia']
    for gene in genes:
        #  Create CDR3 anchors.
        gapped_fasta_file = args.indir + 'abstar_gapped_' + gene + '.fasta'
        gene_dict = get_gene_dict(gapped_fasta_file, gene=gene)
        for software in softwares:
            output_fasta_file_anchor = (args.indir + software + '_'
                                        + gene + '_gene_CDR3_anchors.csv')
            write_to_csv(output_fasta_file_anchor, gene_dict, program=software)
        #  Create genomic references.
        genomic_file = args.indir + 'abstar_genomic_' + gene + 's.fasta'
        get_genomic_ref(genomic_file,
                        genomic_file.replace(".fasta", "_for_igor.fasta"),
                        gene_dict)

if __name__ == '__main__':
    main()
