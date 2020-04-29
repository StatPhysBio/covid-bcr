from utils import *
import csv

def abstar_naive_productive(abstar_file):
    sonia_input = []
    print("Looking at", abstar_file)
    annotations = unpickle(abstar_file)
    for ann in annotations:
        #  Look at naive productive sequences
        if ann['mut_count_nt'] == 0 and ann['prod'] == 'yes':
            #  Skip sequences with N in the middle
            if remove_N_buffer(ann['raw_input']).count("N") != 0:
                continue
            #  Skip sequences with shm indels.
            if ann['vdj_nt'] != ann['gapped_vdj_nt'] or ann['vdj_germ_nt'] != ann['gapped_vdj_germ_nt'] :
                continue
            cdr3 = ann['junc_nt']
            v_gene = ann['v_gene']['gene']
            j_gene = ann['j_gene']['gene']
            time = get_time(ann['seq_id'])
            #abundance = get_dupcounts(ann['seq_id'])
            sonia_input.append((cdr3, v_gene, j_gene, time))

    # Keep only unique nucleotide arrangements
    sonia_input = list(set(sonia_input))
    #  Translate CDR3 from nucleotide to amino acid
    sonia_input = [(translate(si[0]), si[1], si[2], si[3])
                   for si in sonia_input]
    return sonia_input

def partis_naive_productive(partis_file):
    sonia_input = []
    print("Looking at", partis_file)
    annotations = unpickle(partis_file)
    for ann in annotations:
        if ann['invalid']:
            continue
        #  Naive productive sequences.
        if ann['n_mutations'][0] == 0:
            #  Skip if has shm indels, stop codons, or is out of frame.
            if ann['has_shm_indels'][0] or ann['stops'][0] or not ann['in_frames'][0]:
                continue
            #  Skip sequences with N in the middle.
            if remove_N_buffer(ann['input_seqs'][0]).count("N") != 0:
                continue
            cdr3_nt = ann['cdr3_seqs'][0]
            v_gene = ann['v_gene'].split("*")[0]
            j_gene = ann['j_gene'].split("*")[0]
            time = get_time(ann['unique_ids'][0])
            #abundance = get_dupcounts(ann['unique_ids'][0])
            sonia_input.append((cdr3_nt, v_gene, j_gene, time))

    # Keep only unique nucleotide arrangements
    sonia_input = list(set(sonia_input))
    #  Translate CDR3 from nucleotide to amino acid
    sonia_input = [(translate(si[0]), si[1], si[2], si[3])
                   for si in sonia_input]
    return sonia_input

def save_to_csv(outfile, sonia_input):
    if ".csv" not in outfile:
        outfile += ".csv"
    col_titles = 'amino_acid,v_gene,j_gene,time\n'
    with open(outfile, 'w') as f:
        f.write(col_titles)
        for s_in in sonia_input:
            f.write("%s,%s,%s,%s\n"%(s_in[0], s_in[1], s_in[2], s_in[3]))
    print("Wrote sonia input to", outfile)

def main():
    import argparse
    parser = argparse.ArgumentParser(
        description='Get input for SONIA from pickled annotations')
    parser.add_argument('--infile', type=str, help='path to pickled annotations')
    parser.add_argument('--partis', action='store_true', help='annotations from partis')
    parser.add_argument('--abstar', action='store_true', help='annotations from abstar')
    parser.add_argument('--outfile', type=str, help='path to save csv')

    args = parser.parse_args()
    if args.partis:
        sonia_input = partis_naive_productive(args.infile)
    if args.abstar:
        sonia_input = abstar_naive_productive(args.infile)

    save_to_csv(args.outfile, sonia_input)

if __name__ == '__main__':
    main()
