#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thur May 28 2020

@author: Zachary Montague
"""
fasta_file = '/gscratch/stf/zachmon/covid/covid-bcr/abstar_gapped_J.fasta' #name of IMGT generated fasta file with gaps

output_fasta_file_anchor = 'igor_J_gene_CDR3_anchors.csv' #name of output file
with open(fasta_file) as infile:
        with open(output_fasta_file_anchor, 'w') as outfile:
            outfile.write('gene;anchor_index' + '\n')
            for line in infile:
                if line[0] == '>':
                    header = line.split("|")
                    s = header[12] #column with the gene length and number of gaps
                    gaps = int(s.split("+")[0]) #get number of gaps
                    if header[1][-2:]=='01' or header[1][-2:]=='02': #do this only for the first allele. this can be changed if all alleles are needed
                        outfile.write(header[1] + ';' + str(gaps - 34) + '\n') #write name of gene and anchor position (Cys)
output_fasta_file_anchor = 'sonia_J_gene_CDR3_anchors.csv' #name of output file
with open(fasta_file) as infile:
        with open(output_fasta_file_anchor, 'w') as outfile:
            outfile.write('gene,anchor_index,function' + '\n')
            for line in infile:
                if line[0] == '>':
                    header = line.split("|")
                    s = header[12] #column with the gene length and number of gaps
                    gaps = int(s.split("+")[0]) #get number of gaps
                    if header[1][-2:]=='01' or header[1][-2:]=='02': #do this only for the first allele. this can be changed if all alleles are needed
                        print(gaps)
                        outfile.write(header[1] + ',' + str(gaps - 34) +  ',' + header[3] + '\n') #write name of gene and anchor position (Cys)

