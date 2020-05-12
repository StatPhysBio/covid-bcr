from utils import *

def create_new_fasta(save_name, headers, sequences):
    with open(save_name, "w") as new_fasta:
        for i, header in enumerate(headers):
            new_fasta.write(">"+header + "\n")
            new_fasta.write(sequences[i].replace("N","-") + "\n")

def get_consensus_sequence(lineage):
    naives = [ann['naive_seq'] for ann in lineage]
    naive_lengths = list(set([len(n) for n in naives]))
    most_common_length = Counter([len(n) for n in naives]).most_common()[0][0]

    #  Make all naive seqs the same length
    if len(naive_lengths) > 1:
        #  Get an example for the padding
        for seq in naives:
            if len(seq) == most_common_length:
                most_common_length_naive_seq = seq
                break
        #  Padding for max length
        mc_num_n_start, mc_num_n_end = get_num_n_at_start_and_end(most_common_length_naive_seq)

        for idx, seq in enumerate(naives):
            if len(seq) == most_common_length:
                continue
            naives[idx] = fix_N_padding(seq, mc_num_n_start, mc_num_n_end)

    #  Initialize dictionary to hold residues at each position
    consensus_residue_list = {}
    for i in range(most_common_length):
        consensus_residue_list[i] = []

    #  Get residues at each position
    for seq in naives:
        residue_list = list(seq)
        for i in range(most_common_length):
            consensus_residue_list[i].append(residue_list[i])

    #  Create counters at each position
    consensus_counter = {}
    consensus_sequence = ""
    total = len(naives)
    for i in range(most_common_length):
        consensus_counter[i] = Counter(consensus_residue_list[i])
        #  Normalize so easier to see
        for residue in consensus_counter[i]:
            consensus_counter[i][residue] /= total
        consensus_sequence += consensus_counter[i].most_common(1)[0][0]
    return consensus_sequence

def make_div_by_3(ann):
    cdr3_start = ann['codon_positions']['v']
    cdr3_length = ann['cdr3_length']
    seq_length = len(ann['naive_seq'])
    cdr3_end = cdr3_start + cdr3_length
    len_before_cdr3 = cdr3_start
    len_after_cdr3 = seq_length - cdr3_end
    adjusted_beginning_index = len_before_cdr3 % 3
    adjusted_end_index = seq_length - len_after_cdr3 % 3
    return adjusted_beginning_index, adjusted_end_index

def make_fastas_for_trees(lineage, ranking, patient, partis, abstar, outdir):
    consensus_header = ["s0_0"]
    consensus_sequence = get_consensus_sequence(lineage)

    #  Get indices of sequence substring which has substrings
    #  before and after cdr3 divisible by 3 (makes tree building easier)
    #  Crutical to find a lineage with same length as consensus_sequence
    for ann in lineage:
        if len(ann['naive_seq']) == len(consensus_sequence):
            start_index, end_index = make_div_by_3(ann)
            #  To be used to fix padding of input seqs.
            target_length_input_seq = ann['input_seq']
            target_length = len(target_length_input_seq)
            break

    consensus_seq = [consensus_sequence[start_index:end_index]]

    save_name = "patient-" + patient + "_lineage-" + str(ranking)

    consensus_save_name = save_name + "_consensus_sequence.fasta"
    create_new_fasta(join(outdir, consensus_save_name), consensus_header, consensus_seq)

    if partis:
        #  Check that input sequences are all the same length as consensus sequence
        input_seqs = [ann['input_seqs'][0] for ann in lineage]
        input_seq_lengths = list(set([len(seq) for seq in input_seqs]))
        if len(input_seq_lengths) != 1:
            #  Get padding of input seq that was retrieved earlier
            #  and is same length as consensus seq
            target_num_n_start, target_num_n_end = get_num_n_at_start_and_end(target_length_input_seq)
            #  Correct sequences with padding of differing length
            for idx, seq in enumerate(input_seqs):
                if len(seq) == target_length:
                    continue
                input_seqs[idx] = fix_N_padding(seq, target_num_n_start, target_num_n_end)

        heads = [ann['unique_ids'][0] for ann in lineage]
        #  Index the sequence so substrings before and after cdr3 are divisible by 3
        seqs = [ann['input_seqs'][0][start_index:end_index] for ann in lineage]

    # TODO Finish abstar implementation.
    if abstar:
        heads = [ann['seq_id'] for ann in lineage]
        seqs = [ann['raw_input'] for ann in lineage]
    fasta_save_name = save_name + ".fasta"
    create_new_fasta(join(outdir, fasta_save_name), heads, seqs)

def main():
    import argparse
    parser = argparse.ArgumentParser(
        description='Convert lineages into fasta and get consensus sequence for trees')
    parser.add_argument('--lineages', type=str, help='path to lineages')
    parser.add_argument('--outdir', type=str, help='path to outfile directory')
    parser.add_argument('--partis', action='store_true', help='annotations were from partis')
    parser.add_argument('--abstar', action='store_true', help='annotations were from abstar')
    parser.add_argument('--patient', type=str, help='patient number')
    args = parser.parse_args()
    lineages = unpickle(args.lineages)

    #  Sort lineages by size so file name gives some information
    vjls_sorted = sorted(lineages, key=lambda e: len(lineages[e]),reverse=True)
    for i, vjl in enumerate(vjls_sorted):
        if len(lineages['productive'][vjl]) < 20:
            continue
        make_fastas_for_trees(lineages[vjl], i, args.patient, args.partis, args.abstar, args.outdir)

if __name__ == '__main__':
    main()
