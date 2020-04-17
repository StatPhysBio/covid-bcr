#!/usr/bin/env python

from Bio.SeqIO import parse
import pandas as pd
import numpy as np
from typing import Tuple, TextIO
from tinyalign import hamming_distance
import sys


def fasta_parse(fasta: str):
    """parse FASTA file and return pandas DataFrame. Assumes C primer, V
    primer, and abundance are the 2nd, 3rd, and 4th id fields delimited by "|"

    fasta: path to FASTA file
    """
    fasta_dat = []
    for seq in parse(fasta, 'fasta'):
        id, cprimer, vprimer, abundance = (x.split('=')[-1]
                                           for x in seq.id.split('|'))
        fasta_dat.append([id, str(seq.seq), cprimer, vprimer, int(abundance)])
    df = pd.DataFrame(fasta_dat,
                      columns=('id', 'sequence', 'C primer', 'V primer',
                               'abundance'))
    df['length'] = df.sequence.str.len()
    return df


def df2fasta(df: pd.DataFrame, file: TextIO = sys.stdout):
    """print FASTA to file (default stdout) from DataFrame as parsed by
    fasta_parse. FASTA sequence is NOT wrapped to 80 characters.
    """
    try:
        for idx in df.index:
            print(f'>{df.loc[idx, "id"]}|CPRIMER={df.loc[idx, "C primer"]}'
                  f'|VPRIMER={df.loc[idx, "V primer"]}'
                  f'|DUPCOUNT={df.loc[idx, "abundance"]}', file=file)
            print(df.loc[idx, 'sequence'], file=file)
    except BrokenPipeError:
        pass


def df2parents(df: pd.DataFrame, file: TextIO = sys.stdout):
    """print two columns: child and parent
    """
    try:
        for idx in df.index:
            for child in df.loc[idx, 'children']:
                print(f'{child}\t{df.loc[idx, "id"]}', file=file)
    except BrokenPipeError:
        pass


def error_correct(df: pd.DataFrame,
                  delta_r: np.float32 = 1,
                  delta_a: np.float32 = 1) -> Tuple[np.ndarray, np.ndarray]:
    """correct errors by clustering with sparse distance computations

    df: DataFrame as returned by fasta_parse
    delta_r: marginal Hamming distance tolerance per decade in log ratio
             abundances
    delta_a: marginal abundance tolerance of clusterable sequences per decade
             in log ratio abundances (default 1)
    """
    df = df.sort_values(by=['V primer', 'C primer', 'length', 'abundance'],
                        ascending=(True, True, False,
                                   False)).reset_index(drop=True)
    assert len(df['C primer'].unique()) == 1

    df['children'] = [[] for _ in range(len(df))]

    parent_idxs = np.where(df.abundance >= 10 ** (1 / delta_r))[0]

    n_clustered = 0

    for ct, i in enumerate(parent_idxs, 1):
        if df.abundance.values[i] == 0:
            continue
        block_idxs = np.where(np.logical_and(
                              df.length.values == df.length.values[i],
                              df['V primer'].values
                              == df['V primer'].values[i]))[0]
        for j in reversed(block_idxs):
            if j == i:
                break
            if (df.abundance.values[i] == df.abundance.values[j]
                    or df.abundance.values[j] == 0):
                continue
            abundance_log_ratio = (np.log10(df.abundance.values[i])
                                   - np.log10(df.abundance.values[j]))
            if 1 / abundance_log_ratio > delta_r:
                break
            if df.abundance.values[j] / abundance_log_ratio > delta_a:
                break
            d = hamming_distance(df.sequence.values[i], df.sequence.values[j])
            if d / abundance_log_ratio <= delta_r:
                df.abundance.values[i] += df.abundance.values[j]
                df.children.values[i].append(df.id.values[j])
                df.abundance.values[j] = 0
                n_clustered += 1
        print(f'{ct / len(parent_idxs):.2%}, corrected {n_clustered}',
              end='    \r', flush=True, file=sys.stderr)
    print(file=sys.stderr)

    return df[df.abundance > 0]


def main():
    """
    usage: python error_correct.py -h"""
    import argparse
    import pickle

    parser = argparse.ArgumentParser(
        description='FASTA error correction, streams to stdout')
    parser.add_argument('fasta', type=str, help='path to FASTA')
    parser.add_argument('outbase', type=str, help='basename for output files')
    parser.add_argument('--delta_r', type=float, default=1,
                        help='marginal Hamming distance tolerance per decade '
                             'in log ratio abundances (default 1)')
    parser.add_argument('--delta_a', type=float, default=1,
                        help='marginal abundance tolerance of clusterable '
                             'sequences per decade in log ratio abundances '
                             '(default 1)')
    parser.add_argument('--passes', type=int, default=1,
                        help='number of times to repeat greedy clustering '
                             '(default 1)')
    parser.add_argument('--keep_singletons', action='store_true',
                        help="don't discard uncorrected singletons")
    args = parser.parse_args()

    df = fasta_parse(args.fasta)

    print("Args chosen:\nfasta",args.fasta,"\nargs.outbase",args.outbase,"\ndel_r", args.delta_r, "\ndel_a", args.delta_a)

    for this_pass in range(1, args.passes + 1):
        if args.passes > 1:
            print(f'pass {this_pass}', file=sys.stderr)
        df = error_correct(df, delta_r=args.delta_r, delta_a=args.delta_a)

    if not args.keep_singletons:
        df = df[df.abundance > 1]

    df2fasta(df, file=open(f'{args.outbase}.corrected.fa', 'w'))
    df2parents(df, file=open(f'{args.outbase}.parents.tsv', 'w'))
    pickle.dump(df, open(f'{args.outbase}.df.pkl', 'wb'))

if __name__ == '__main__':
    main()
