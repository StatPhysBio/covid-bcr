import sonia
from sonia.sonia_vjl import SoniaVJL
from sonia.sonia_leftpos_rightpos import SoniaLeftposRightpos
from sonia.sonia_length_pos import SoniaLengthPos
import pandas as pd

def main():
    import argparse
    parser = argparse.ArgumentParser(
        description='create sonia model')
    parser.add_argument('--sonia_dir', type=str,
                        help='path to sonia model for saving')
    parser.add_argument('--sonia_model', type=str,
                        help='options: [vjl, length, leftright]\n'
                        'default: vjl')
    parser.add_argument('--include_joint_genes', action='store_true',
                        help='Have features for combinations of vj\n'
                        'otherwise independent features for v and j.\n'
                        'default: false')
    parser.add_argument('--in_data', type=str,
                       help='path to data csv file')
    parser.add_argument('--in_gen', type=str,
                        help='path to gen csv file')
    parser.add_argument('--lineage_size', type=int,
                        help='exclusive cutoff for minimum size of lineage')
    parser.add_argument('--epochs', type=int,
                        help='number of epochs. Default=50')
    parser.add_argument('--num_gen', type=int,
                        help='number of gen\n'
                        'Default to generate is 2e5.\n'
                        'Default of in_gen is all seqs.\n')
    args = parser.parse_args()

    models = ['vjl', 'length', 'leftright']
    print(args.sonia_model)
    if args.sonia_model in models:
        model = args.sonia_model
        print(model)
    else:
        model = 'vjl'
    if args.include_joint_genes is not None:
        include_joint_genes = args.include_joint_genes
        include_indep_genes = not args.include_joint_genes
    else:
        include_join_genes = False
        include_indep_genes = True
    if args.lineage_size is None:
        min_size = 0
    else:
        min_size = args.lineage_size
    print("joint genes", include_joint_genes)
    print("indep genes", include_indep_genes)
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
        if model is 'vjl':
            qm = SoniaVJL(data_seqs=data,gen_seqs=gen[:num_gen_seqs],chain_type='humanIGH',
                         include_indep_genes=include_indep_genes, include_joint_genes=include_joint_genes)
        elif model is 'length':
            qm = SoniaLengthPos(data_seqs=data,gen_seqs=gen[:num_gen_seqs],chain_type='humanIGH',
                         include_indep_genes=include_indep_genes, include_joint_genes=include_joint_genes)
        elif model is 'leftright':
            qm = SoniaLeftposRightPos(data_seqs=data,gen_seqs=gen[:num_gen_seqs],chain_type='humanIGH',
                     include_indep_genes=include_indep_genes, include_joint_genes=include_joint_genes)
    else:
        if model == 'vjl':
            qm = SoniaVJL(data_seqs=data, chain_type='humanIGH',
                         include_indep_genes=include_indep_genes, include_joint_genes=include_joint_genes)
        elif model == 'length':
            qm = SoniaLengthPos(data_seqs=data, chain_type='humanIGH',
                         include_indep_genes=include_indep_genes, include_joint_genes=include_joint_genes)
        elif model == 'leftright':
            qm = SoniaLeftposRightpos(data_seqs=data, chain_type='humanIGH',
                         include_indep_genes=include_indep_genes, include_joint_genes =include_joint_genes)
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
