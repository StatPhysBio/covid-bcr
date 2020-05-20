import os
import sonia
from sonia.sonia_vjl import SoniaVJL
from sonia.evaluate_model import EvaluateModel
from sonia.sequence_generation import SequenceGeneration

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

def save_output(save_name, out_seqs, Q=None, pgen=None, ppost=None):
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
    parse.add_argument('--outname', type=str,
                       help='path to output file')
    parser.add_argument('--num_seqs', type=int,
                        help='number of sequences to generate')
    parser.add_argument('--evaluate', action='store_true',
                        help='evaluate sequences to get pgen, ppost, and Q')
    parser.add_argument('--post', action='store_true',
                        help='generate post sequences')
    parser.add_argument('--pre', action='store_true',
                        help='generate pre sequences')

    args = parser.parse_args()
    save_dir = args.sonia_dir

    qm = SoniaVJL(load_dir=args.sonia_dir, chain_type='humanIGH')
    gn = SequenceGeneration(qm)

    if args.post:
        post_seqs = generate_post_seqs(gn, numseqs=args.num_seqs)
    if args.pre:
        pre_seqs = generate_pre_seqs(gn, numseqs=args.num_seqs)

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
    else:
        if args.post:
            save_name = os.path.join(save_dir, 'post_seqs.csv')
            save_output(save_name, post_seqs)
        if args.pre:
            save_name = os.path.join(save_dir, 'pre_seqs.csv')
            save_output(save_name, pre_seqs)
if __name__=='__main__':
    main()
