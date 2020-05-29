"""
Created on Tue Sep 11 2018

@author: Yuval Elhanati
Updated May 28 2020, Zachary Montague
"""

from collections import OrderedDict
from operator import itemgetter
from Bio.SeqIO import parse

def sort_dict(unsorted_dict):
    D = OrderedDict(sorted(unsorted_dict.items(), key=itemgetter(0)))
    return D

def get_v_gaps(header):
    col = header[12]
    return int(col[col.find('+')+1:col.find('=')])

def get_v_anchor(gaps):
    return str(3*(104-1) - gaps)

def get_j_gaps(header):
    col = header[12]
    gaps = int(col.split("+")[0])
    return gaps

def get_j_anchor(gaps):
    return str(gaps - 34)

def trim_dict(gene_dict):
    gene_dict = sort_dict(gene_dict)
    out_dict = {}
    for gene in gene_dict:
        gene_dict[gene] = sort_dict(gene_dict[gene])
        for idx, allele in enumerate(gene_dict[gene]):
            if idx > 1:
                break
            out_dict[gene + "*" + allele] = gene_dict[gene][allele]
    return out_dict

def get_gene_dict(gapped_fasta, gene='V'):
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
                line = line.split("|")
                gene = line[1].split("*")[0]
                if gene not in gene_dict:
                    gene_dict[gene] = {}
                allele = line[1].split("*")[-1]
                gaps = gap_func(line)
                gene_dict[gene][allele] = anchor_func(gaps)
    return trim_dict(gene_dict)

def write_to_csv(save_name, gene_dict, delimiter):
    with open(save_name, 'w') as outfile:
        outfile.write('gene' + delimiter + 'anchor_index\n')
        for key in gene_dict:
            outfile.write(key + delimiter + gene_dict[key] + "\n")

def get_genomic_ref(fasta, save_name, gene_dict):
    seqs = []
    headers = []
    for seq in parse(fasta, 'fasta'):
        headers.append(seq.description)
        seqs.append(str(seq.seq))

    with open(save_name, 'w') as f:
        for idx,header in enumerate(headers):
            if header in gene_dict:
                f.write(">"+header+"\n")
                f.write(seqs[idx] + "\n")

def main():
    genes = ['V', 'J']
    out_dict = {'igor': ';', 'sonia':','}
    for gene in genes:
        gapped_fasta_file = '/gscratch/stf/zachmon/covid/covid-bcr/abstar_gapped_'+gene+'.fasta'
        gene_dict = get_gene_dict(gapped_fasta_file, gene=gene)
        for key in out_dict:
            output_fasta_file_anchor = key + '_' + gene + '_gene_CDR3_anchors.csv'
            write_to_csv(output_fasta_file_anchor, gene_dict, out_dict[key])
        genomic_file = '/gscratch/stf/zachmon/covid/covid-bcr/abstar_genomic_' + gene + 's.fasta'
        get_genomic_ref(genomic_file,
                        genomic_file.replace(".fasta", "_for_igor.fasta"),
                        gene_dict)

if __name__ == '__main__':
    main()
