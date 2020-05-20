import sonia
from sonia.sonia_vjl import SoniaVJL
import pandas as pd

def main():
    import argparse
    parser = argparse.ArgumentParser(
        description='create sonia model')
    parser.add_argument('--sonia_dir', type=str,
                        help='path to sonia model for saving')
    parser.add_argument('--in_data', type=str,
                       help='path to data csv file')
    parser.add_argument('--in_gen', type=str,
                        help='path to gen csv file')
    parser.add_argument('--lineage_size', type=int,
                        help='exclusive cutoff for minimum size of lineage')
    parser.add_argument('--epochs', type=int,
                        help='number of epochs. Default=50')
    parser.add_argument('--num_gen', type=int,
                        help='number of gen. Default to generate is 2e5. Default of in_gen is all seqs.')
    args = parser.parse_args()
    if args.lineage_size is None:
        min_size = 0
    else:
        min_size = args.lineage_size
    data_df = pd.read_csv(args.in_data)
    trimmed_df = data_df[data_df['lineage_size'] > min_size]
    data = trimmed_df[['consensus_cdr3','v_gene','j_gene']].values.tolist()
    print("Amount of data:", len(data))

    #  Upper bound for the amount of gen seqs.
    #  Otherwise the probability of not having
    #  a data seq in a batch (size 5e3 by default)
    #  is non-neglible.
    num_gen_seqs = 100*len(data)

    if args.in_gen is not None:
        gen = pd.read_csv(args.in_gen)[['amino_acid','v_gene','j_gene']].values.tolist()
        if args.num_gen is not None:
            if args.num_gen < num_gen_seqs:
                num_gen_seqs = args.num_gen
        print("Amount of gen:",num_gen_seqs)
        qm = SoniaVJL(data_seqs=data,gen_seqs=gen[:num_gen_seqs],chain_type='humanIGH',
                     include_indep_genes = True, include_joint_genes = False )
    else:
        qm = SoniaVJL(data_seqs=data, chain_type='humanIGH',
                     include_indep_genes = True, include_joint_genes = False )
        if args.num_gen is not None:
            if args.num_gen < num_gen_seqs:
                num_gen_seqs = args.num_gen
            print("Amount of gen:",num_gen_seqs)
            qm.add_generated_seqs(args.num_gen)
        else:
            num_gen_seqs = int(2e5)
            print("Amount of gen:",num_gen_seqs)
            qm.add_generated_seqs(num_gen_seqs)
    if args.epochs is None:
        epochs = 50
    else:
        epochs = args.epochs
    qm.infer_selection(epochs=epochs)
    qm.save_model(args.sonia_dir)

if __name__=='__main__':
    main()
