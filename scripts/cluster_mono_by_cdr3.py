from abstar_pipeline import *
from vjl_slc import *
from utils import *
from jellyfish import hamming_distance
import pandas as pd
import numpy as np

def make_v_ref_dict():
    v_anchors = pd.read_csv('/gscratch/stf/zachmon/covid/covid-bcr/sonia_input/V_gene_CDR3_anchors.csv')
    v_anchors=v_anchors.assign(g=v_anchors['gene'].str.split("*").str[0])
    anchor_ind = v_anchors.drop_duplicates('g').set_index('g').loc[v_anchors['g'].unique()]['anchor_index']
    vs,heads = fasta_read('/gscratch/stf/zachmon/covid/covid-bcr/igor_input/abstar_genomic_Vs.fasta')
    v_dict = {}
    for idx,h in enumerate(heads):
        gene = h.split("*")[0]
        if gene in v_dict:
            continue
        v_trim = vs[idx][:anchor_ind.loc[gene]]
        if v_trim[-3:] == 'TGT':
            v_trim = v_trim[:-3]
        v_dict[gene] = v_trim
    v_dict = sort_dict(v_dict)
    return v_dict

v_ref = make_v_ref_dict()
v_ref_genes = list(v_ref.keys())
v_ham_mat = np.zeros(shape=(len(v_ref),len(v_ref)))

for i,v1 in enumerate(v_ref):
    for j,v2 in enumerate(v_ref):
        g1 = v_ref[v1]
        g2 = v_ref[v2]
        min_len = np.min([len(g1),len(g2)])
        v_ham_mat[i,j] = hamming_distance(g1[-min_len:],g2[-min_len:]) / min_len

v_to_include = {}
for i,v1 in enumerate(v_ref_genes):
    v_to_include[v1] = []
    for j,v2 in enumerate(v_ref_genes):
        if v_ham_mat[i,j] <= 0.1:
            v_to_include[v1].append(v2)


def cluster_by_l(ann, den, lins, thresh):
    v_gene = ann['v_gene']['gene']
    cdr3 = ann['junc_nt']
    len_cdr3 = len(cdr3)
    subset = [key for key in den
              if key[2] == str(len_cdr3)
              and key[0] in v_to_include[v_gene]]
    test = []
    for j,k in enumerate(subset):
        testk = np.zeros(len(den[k]),dtype=np.float16)
        for i,a in enumerate(den[k]):
            testk[i] = hamming_distance(cdr3, a['junc_nt'])
        testk /= len_cdr3
        test.append(np.min(testk))
    test = np.array(test)
    to_insert = np.array(subset)[test < thresh]
    for key in to_insert:
        print(ann['seq_id'].split('|')[0], ann['seq_id'].split('|')[-1], key)
        lins['productive'][key[0]][key[1]][key[2]][key[3]] += [ann]

def cluster_by_l_aa(ann, den, lins, thresh):
    v_gene = ann['v_gene']['gene']
    cdr3 = ann['junc_nt']
    cdr3_aa = translate(cdr3)
    len_cdr3 = len(cdr3)
    subset = [key for key in den
              if key[2] == str(len_cdr3)
              and key[0] in v_to_include[v_gene]]
    test = []
    for j,k in enumerate(subset):
        testk = np.zeros(len(den[k]),dtype=np.float16)
        for i,a in enumerate(den[k]):
            testk[i] = hamming_distance(cdr3_aa, translate(a['junc_nt']))
        testk /= len(cdr3_aa)
        test.append(np.min(testk))
    test = np.array(test)
    to_insert = np.array(subset)[test < thresh]
    for key in to_insert:
        print(ann['seq_id'].split('|')[-1], key)
        lins['productive'][key[0]][key[1]][key[2]][key[3]] += [ann]

def main():
    import argparse
    parser = argparse.ArgumentParser(
        description='')
    parser.add_argument('--infile', type=str,
                        help='.json lineage files. e.g. usage: /PATH/TO/LINEAGES/*')
    args = parser.parse_args()

    mono = json_open('/gscratch/stf/zachmon/covid/monoclonal_ntd_rbd.json')
    lins = json_open(args.infile)
    den = denest_lineages(lins['productive'])

    for ann in mono:
        cluster_by_l(ann, den, lins, 0.15)
    #save_dir = '/gscratch/spe/zacmon/bulk_plasma_rbd_ntd/cdr3_nt_cluster/'
    #save_dir = '/gscratch/stf/zachmon/covid/11_3_lineages_bulk_plasma_rbd_ntd/sc_lins/'
    save_dir = '/gscratch/stf/zachmon/covid/bullshit/'
    save_name = args.infile.split('/')[-1].replace('.json','_with_sc.json')
    json_save(save_dir + save_name, lins)

if __name__ == '__main__':
    main()
