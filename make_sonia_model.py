import sonia
from sonia.sonia_vjl import SoniaVJL
import pickle
import pandas as pd


def get_patient(uid):
    return int(uid.split("|")[0].split("-")[-1])

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
                        help='cutoff for minimum size of lineage')
    parser.add_argument('--epochs', type=int,
                        help='number of epochs. Default=50')
    parser.add_argument('--num_gen', type=int,
                        help='number of gen. Default=2e6')
    args = parser.parse_args()
    gen = pd.read_csv(args.in_gen)[['amino_acid','v_gene','j_gene']].values.tolist()
    if args.num_gen is not None:
        gen = gen[:args.num_gen]
    if args.lineage_size is None:
        min_size = 0
    else:
        min_size = args.lineage_size
    data_df = pd.read_csv(args.in_data)
    trimmed_df = data_df[data_df['lineage_size'] > min_size]
    data = trimmed_df[['consensus_cdr3','v_gene','j_gene']].values.tolist()

    qm = SoniaVJL(data_seqs=data,gen_seqs=gen,chain_type='humanIGH')
    if args.epochs is None:
        epochs = 50
    else:
        epochs = args.epochs
    qm.infer_selection(epochs=epochs)
    qm.save_model(args.sonia_dir)

if __name__=='__main__':
    main()
