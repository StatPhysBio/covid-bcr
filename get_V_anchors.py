#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 11 2018

@author: Yuval Elhanati
"""
fasta_file = '/Users/zach/Documents/github/covid-bcr/abstar_gapped_IGH.fasta' #name of IMGT generated fasta file with gaps

output_fasta_file_anchor = 'V_gene_CDR3_anchors.csv' #name of output file
with open(fasta_file) as infile:
        with open(output_fasta_file_anchor, 'w') as outfile:
            outfile.write('gene,anchor_index' + '\n')
            for line in infile:
                if line[0] == '>':
                    header = line.split("|")
                    s = header[12] #column with the gene length and number of gaps
                    gaps = int(s[s.find('+')+1:s.find('=')]) #get number of gaps
                    if header[1][-2:]=='01': #do this only for the first allele. this can be changed if all alleles are needed
                        outfile.write(header[1] + ',' + str(3*(104-1) - gaps) + '\n') #write name of gene and anchor position (Cys)
