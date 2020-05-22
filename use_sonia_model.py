import os
import sonia
from sonia.sonia_vjl import SoniaVJL
from sonia.evaluate_model import EvaluateModel
from sonia.sequence_generation import SequenceGeneration
import pandas as pd
def generate_pre_seqs(gn, numseqs=None):
    if numseqs is None:
        numseqs=10000
    pre_seqs=gn.generate_sequences_pre(numseqs)
    return pre_seqs

def generate_post_seqs(gn, numseqs=None):
    if numseqs is None:
        numseqs=10000
    post_seqs=gn.generate_sequences_post(numseqs)
    return post_seqs

def evaluate_seq(ev, in_seqs):
    Q_model,pgen_model,ppost_model=ev.evaluate_seqs(in_seqs)
    return Q_model, pgen_model, ppost_model

def save_output(save_name, out_seqs, data=False, Q=None, pgen=None, ppost=None):
    if data:
        with open(save_name,'w') as seqs_file:
            seqs_file.write('cdr3,vgene,jgene,lineage_size,expanded,Q,pgen,ppost\n')
            for i,seq in enumerate(out_seqs):
                seqs_file.write("%s,%s,%s,%i,%s,%.16e,%.16e,%.16e\n"%(seq[0],seq[1],seq[2],seq[3],seq[4],Q[i],pgen[i],ppost[i]))
        print("Wrote",save_name)
        return
    if Q is None:
        with open(save_name, 'w') as seqs_file:
            seqs_file.write('cdr3,vgene,jgene\n')
            for i,seq in enumerate(out_seqs):
                seqs_file.write("%s,%s,%s\n"%(seq[0],seq[1],seq[2]))
        print("Wrote",save_name)
    else:
        with open(save_name,'w') as seqs_file:
            seqs_file.write('cdr3,vgene,jgene,Q,pgen,ppost\n')
            for i,seq in enumerate(out_seqs):
                seqs_file.write("%s,%s,%s,%.16e,%.16e,%.16e\n"%(seq[0],seq[1],seq[2],Q[i],pgen[i],ppost[i]))
        print("Wrote",save_name)

def main():
    import argparse
    parser = argparse.ArgumentParser(
        description='use sonia model to generate and evaluate sequences')
    parser.add_argument('--sonia_dir', type=str,
                        help='path to sonia model')
    parser.add_argument('--out_dir', type=str,
                       help='path to output')
    parser.add_argument('--num_seqs', type=int,
                        help='number of sequences to generate')
    parser.add_argument('--in_data', type=str,
                        help='input file of csv of seqs, v, and j')
    parser.add_argument('--evaluate', action='store_true',
                        help='evaluate sequences to get pgen, ppost, and Q')
    parser.add_argument('--post', action='store_true',
                        help='generate post sequences')
    parser.add_argument('--pre', action='store_true',
                        help='generate pre sequences')
    parser.add_argument('--lineage_size', type=int,
                        help='exclusive cutoff on lineage size')
    parser.add_argument('--evaluate_data_for_model', action='store_true',
                        help='evaluate data used to generate model')
    parser.add_argument('--evaluate_gen_for_model', action='store_true',
                        help='evaluate gen used to generate model')
    args = parser.parse_args()
    if args.lineage_size is None:
        min_size = 0
    else:
        min_size = args.lineage_size
    if args.in_data is not None:
        data_df = pd.read_csv(args.in_data)
        trimmed_df = data_df[data_df['lineage_size'] > min_size]
        all_data = trimmed_df[['consensus_cdr3', 'v_gene', 'j_gene', 'lineage_size', 'expanded']].values.tolist()
        data = trimmed_df[['consensus_cdr3','v_gene','j_gene']].values.tolist()


    save_dir = args.sonia_dir
    if args.out_dir is not None:
        save_dir = args.out_dir
    qm = SoniaVJL(load_dir=args.sonia_dir, chain_type='humanIGH')
    gn = SequenceGeneration(qm)

    if args.post:
        post_seqs = generate_post_seqs(gn, numseqs=args.num_seqs)
    if args.pre:
        pre_seqs = generate_pre_seqs(gn, numseqs=args.num_seqs)

    if args.evaluate_data_for_model:
        ev = EvaluateModel(qm)
        Q_out, pgen_out, ppost_out = evaluate_seq(ev, qm.data_seqs)
        save_name = os.path.join(save_dir, 'model_data_eval.csv')
        save_output(save_name, qm.data_seqs, Q=Q_out, pgen=pgen_out, ppost=ppost_out)
    if args.evaluate_gen_for_model:
        ev = EvaluateModel(qm)
        if args.num_seqs is not None:
            num_to_eval = args.num_seqs
        else:
            num_to_eval = len(qm.gen_seqs)
        Q_out, pgen_out, ppost_out = evaluate_seq(ev, qm.gen_seqs[:num_to_eval])
        save_name = os.path.join(save_dir, 'model_gen_eval.csv')
        save_output(save_name, qm.gen_seqs[:num_to_eval], Q=Q_out, pgen=pgen_out, ppost=ppost_out)

    if args.evaluate:
        ev = EvaluateModel(qm)
        if args.post:
            Q_out, pgen_out, ppost_out = evaluate_seq(ev, post_seqs)
            save_name = os.path.join(save_dir, 'post_seqs_and_eval.csv')
            save_output(save_name, post_seqs, Q=Q_out, pgen=pgen_out, ppost=ppost_out)
        if args.pre:
            Q_out, pgen_out, ppost_out = evaluate_seq(ev, pre_seqs)
            save_name = os.path.join(save_dir, 'pre_seqs_and_eval.csv')
            save_output(save_name, pre_seqs, Q=Q_out, pgen=pgen_out, ppost=ppost_out)
        if args.in_data is not None:
            Q_out, pgen_out, ppost_out = evaluate_seq(ev, data)
            save_name = os.path.join(save_dir, 'data_eval.csv')
            save_output(save_name, all_data, data=True,Q=Q_out, pgen=pgen_out, ppost=ppost_out)
    else:
        if args.post:
            save_name = os.path.join(save_dir, 'post_seqs.csv')
            save_output(save_name, post_seqs)
        if args.pre:
            save_name = os.path.join(save_dir, 'pre_seqs.csv')
            save_output(save_name, pre_seqs)
if __name__=='__main__':
    main()
