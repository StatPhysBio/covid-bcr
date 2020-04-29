#!/usr/bin/env python

from Bio.SeqIO import parse
import pandas as pd
import numpy as np
from typing import Tuple, TextIO
from jellyfish import hamming_distance
import multiprocessing


def fasta_parse(fasta: str):
    """parse FASTA file and return pandas DataFrame. Assumes C primer, V
    primer, and abundance are the 2nd, 3rd, and 4th id fields delimited by "|"

    fasta: path to FASTA file
    """
    fasta_dat = []
    index = []
    for seq in parse(fasta, 'fasta'):
        id, cprimer, vprimer, abundance = (x.split('=')[-1]
                                           for x in seq.id.split('|'))
        index.append(id)
        fasta_dat.append([str(seq.seq), cprimer, vprimer, int(abundance)])

    df = pd.DataFrame(fasta_dat,
                      index=index,
                      columns=('sequence', 'C primer', 'V primer',
                               'abundance'))
    df['length'] = df.sequence.str.len()
    return df


def df2fasta(df: pd.DataFrame, file: TextIO):
    """print FASTA to file from DataFrame as parsed by
    fasta_parse. FASTA sequence is NOT wrapped to 80 characters.
    """
    try:
        for idx in df.index:
            print(f'>{idx}|CPRIMER={df.loc[idx, "C primer"]}'
                  f'|VPRIMER={df.loc[idx, "V primer"]}'
                  f'|DUPCOUNT={df.loc[idx, "abundance"]}', file=file)
            print(df.loc[idx, 'sequence'], file=file)
    except BrokenPipeError:
        pass


def df2parents(df: pd.DataFrame, file: TextIO):
    """print two columns: child and parent
    """
    try:
        for idx in df.index:
            for child in df.loc[idx, 'children']:
                print(f'{child}\t{idx}', file=file)
    except BrokenPipeError:
        pass


def error_correct(df: pd.DataFrame, d_tol: int = 1, a_tol: np.float32 = None,
                  name: str = '', verbose: bool = False
                  ) -> Tuple[np.ndarray, np.ndarray]:
    if a_tol is None:
        a_tol = 1

    # sort by descending abundance
    df = df.sort_values(by='abundance', ascending=False)

    # add a field for the list of children (if any), initialized empty
    df['children'] = [[] for _ in range(len(df))]

    # the possible parents have abundance at least a_tol, since minimum
    # abundance of child is 1
    parent_idxs = np.where(df.abundance.values >= a_tol)[0]

    n_clustered = 0

    for ct, i in enumerate(parent_idxs, 1):
        if df.abundance.values[i] == 0:
            continue
        for j in reversed(range(len(df))):
            if j == i:
                break
            if df.abundance.values[j] == 0:
                continue
            abundance_ratio = df.abundance.values[i] / df.abundance.values[j]
            if abundance_ratio < a_tol:
                break
            if hamming_distance(df.sequence.values[i],
                                df.sequence.values[j]) <= d_tol:
                df.abundance.values[i] += df.abundance.values[j]
                df.children.values[i].append(df.index.values[j])
                df.abundance.values[j] = 0
                n_clustered += 1
        if verbose:
            print(name, f'{ct / len(parent_idxs):.2%}', end='      \r')
    if verbose:
        print()

    return df[df.abundance > 0]


def main():
    """
    usage: python error_correct.py -h"""
    import argparse
    import pickle

    parser = argparse.ArgumentParser(
        description='FASTA error correction')
    parser.add_argument('fasta', type=str, help='path to FASTA')
    parser.add_argument('outbase', type=str, help='basename for output files')
    parser.add_argument('--d_tol', type=float, default=1,
                        help='Hamming distance tolerance (default 1)')
    parser.add_argument('--a_tol', type=float, default=None,
                        help='abundance ratio tolerance (if None, not '
                             'applied)')
    parser.add_argument('--jobs', type=int, default=1,
                        help='number of parallel jobs (default 1)')
    parser.add_argument('--passes', type=int, default=1,
                        help='number of times to repeat greedy clustering '
                             '(default 1)')
    parser.add_argument('--keep_singletons', action='store_true',
                        help="don't discard uncorrected singletons")
    parser.add_argument('--verbose', action='store_true',
                        help='print progress messages')
    args = parser.parse_args()

    df = fasta_parse(args.fasta)

    for this_pass in range(1, args.passes + 1):
        if args.passes > 1 and args.verbose:
            print(f'pass {this_pass}')
        groups = df.groupby(by=['V primer', 'C primer', 'length'])
        with multiprocessing.Pool(processes=args.jobs) as pool:
            chunk_gen = ((df_group, args.d_tol, args.a_tol, name, args.verbose)
                         for name, df_group in groups)
            df = pd.concat(pool.starmap(error_correct, chunk_gen))

    if not args.keep_singletons:
        df = df[df.abundance > 1]

    df2fasta(df, open(f'{args.outbase}.corrected.fa', 'w'))
    df2parents(df, open(f'{args.outbase}.parents.tsv', 'w'))
    pickle.dump(df, open(f'{args.outbase}.df.pkl', 'wb'))


if __name__ == '__main__':
    main()
