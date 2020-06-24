from Bio.SeqIO import parse
import numpy as np
from typing import List, Tuple, TextIO
from jellyfish import hamming_distance

def get_cprimer(header: str) -> str:
    """Obtains the cprimer from the header.

    Parameters
    ----------
    header : str

    Returns
    -------
    str
        cprimer
    """

    return header.split("|")[1].split("=")[-1]

def get_vprimer(header: str) -> str:
    """Obtains the vprimer from the header.

    Parameters
    ----------
    header : str

    Returns
    -------
    str
        vprimer
    """

    return header.split("|")[2].split("=")[-1]

def get_abundance(header: str) -> int:
    """Obtains the abundance counts from the header.

    Parameters
    ----------
    header : str

    Returns
    -------
    int
        abundance counts
    """

    return int(header.split("|")[3].split("=")[-1])

def get_time(header: str) -> int:
    """Obtains the time from the header.

    Parameters
    ----------
    header : str

    Returns
    -------
    int
        time
    """

    return int(header.split("|")[4].split("=")[-1])

def get_replicate(header: str) -> int:
    """Obtains the replicate from the header.

    Parameters
    ----------
    header : str

    Returns
    -------
    int
        replicate
    """

    return int(header.split("|")[6].split("=")[-1])

def fasta_parse(fasta: str) -> Tuple[List[str], List[str]]:
    """Reads a fasta file and prases the headers and the sequences.

    Parameters
    ----------
    fasta : str
        Fasta file containing sequence information.

    Returns
    -------
    headers :  list
        List of header information from the fasta file.
    sequences : list
        List of sequences from the fasta file.
    """

    headers = []
    sequences = []
    for sequence in parse(fasta, 'fasta'):
        sequences.append(str(sequence.seq))
        headers.append(sequence.id)
    return headers, sequences

def write_to_fasta(outfile: str, headers: List[str], sequences: List[str]) -> None:
    """Writes headers and sequences to a fasta file quickly.

    Parameters
    ----------
    outfile : str
        Path and name of file which will be written.
    headers : list
        List of headers to write to file.
    sequences : list
        List of sequences to write to file.

    Returns
    -------
    None
    """

    with open(outfile, "w") as outf:
        for i, header in enumerate(headers):
            outf.write(">"+header + "\n" + sequences[i] + "\n")

def update_header(header: str, abundance: int) -> str:
    """Changes abundance count in header.

    Parameters
    ----------
    header : str
        String containing information about a sequence/annotation,
    abundance : int
        Abundance count to replace that currently in header.

    Returns
    -------
    updated_header : str
        String containing updated abundance.
    """

    header_split = header.split("|")
    header_split[3] = "DUPCOUNT=" + str(abundance)
    updated_header = "|".join(header_split)
    return updated_header

def sort_data_by_abundance(headers: List[str], sequences: List[str]) -> Tuple[List[str], List[str]]:
    """Sorts headers and sequences by descending abundance.

    Parameters
    ----------
    headers : list
        List of headers to be updated.
    sequences : list
        List of sequences to be processed.

    Returns
    -------
    sorted_headers : list
        Headers sorted by abundance.
    sorted_sequences : list
        Sequences sorted by header abundance.
    """

    sorted_headers = sorted(headers, key=lambda h: get_abundance(h), reverse=True)
    sorted_sequences = [seq
                        for _,seq in sorted(zip(headers,sequences),
                                            key=lambda pair: get_abundance(pair[0]),
                                                reverse=True)]
    return sorted_headers, sorted_sequences

def error_correct_marginal(headers: List[str], sequences: List[str],
                           delta_r: np.float32 = 1.0,
                           delta_a: np.float32 = 1.0,
                           debug: bool = False) -> Tuple[List[str], List[str]]:
    """Corrects sequencing errors that separated large abundance clones.

    Parameters
    ----------
    headers : list
        List of headers.
    sequences : list
        List of sequences which all have the same length.
    delta_r : numpy.float32, optional
        The marginal Hamming distance tolerance per decade in log ratio
        abundance (each log_10 unit allowing delta_r additional
        sequence differences).
    delta_a : numpy.float32, optional
        The marginal abundance tolerance of clusterable sequences per decade
        in log ratio abundance (each log_10 unit allowing abundance Delta_a
        higher as clusterable)
    debug : bool, optional
        Bool used to control whether or not to return parent_children dictionary.

    Returns
    -------
    updated_headers : list
        List of headers updated with merged abundances.
    updated_sequences : list
        List of sequences corresponding to the aforementioned headers.
    parent_child : dict
        Dictionary containing which headers got absorbed into what other headers.
    """

    sorted_headers, sorted_sequences = sort_data_by_abundance(headers, sequences)
    sorted_abundances = [get_abundance(h) for h in sorted_headers]

    #  Before any error correction, delta_r already places a limit on which sequences
    #  can be merged into.
    parent_indices = np.where(np.array(sorted_abundances) >= 10 ** (1 / delta_r))[0]

    #  Dictionary used to debug. In conjunction with abstar or annotation output,
    #  it can be used to track how sequences of different types of productivity
    #  or with/without SHM indels get merged.
    parent_child = {}

    for i in parent_indices:
        #  The sequence has already been absorbed.
        if sorted_abundances[i] == 0:
            continue
        parent_name = sorted_headers[i].split("|")[0]
        for j,seq2 in reversed(list(enumerate(sorted_sequences))):
            if i == j:
                break
            #  This algorithm is meant to absorb sequences into larger clones.
            if (sorted_abundances[i] == sorted_abundances[j]
                or sorted_abundances[j] == 0):
                continue
            abundance_log_ratio = (np.log10(sorted_abundances[i])
                                   - np.log10(sorted_abundances[j]))
            #  Since a sequence gets absorbed if d <= delta_r * abundance_log_ratio,
            #  if 1 > delta_r * abundance_log_ratio, then d <= delta_r * abundance_log_ratio
            #  will never be satisfied. We are assuming d = 0 sequences are deduplicated.
            if 1 / abundance_log_ratio > delta_r:
                break
            #  sorted_abundance[j] < delta_a * abundance_log_ratio for a sequence
            #  to be merged.
            if sorted_abundances[j] / abundance_log_ratio > delta_a:
                break
            seq1 = sorted_sequences[i]
            d = hamming_distance(seq1, seq2)
            #  d <= delta_r * abundance_log_ratio for a sequence to be merged.
            if d / abundance_log_ratio <= delta_r:
                if parent_name not in parent_child:
                    parent_child[parent_name] = []
                #  Add merged sequences as children to parent sorted_headers[i].
                parent_child[parent_name].append(sorted_headers[j].split("|")[0])
                #  Update the abundances.
                sorted_abundances[i] += sorted_abundances[j]
                sorted_abundances[j] = 0

    #  Get resulting sequences and headers
    updated_sequences = []
    updated_headers = []

    for index, abun in enumerate(sorted_abundances):
        updated_sequences.append(sorted_sequences[index])
        updated_headers.append(update_header(sorted_headers[index], sorted_abundances[index]))
    if debug:
        return updated_headers, updated_sequences, parent_child
    else:
        return updated_headers, updated_sequences

def error_correct_total(headers: List[str], sequences: List[str],
                        d_thresh: int = 1,
                        a_thresh: np.float32 = None,
                        debug: bool = False) -> Tuple[List[str], List[str]]:
    """Targets correction of reverse transcriptase errors.

    Parameters
    ----------
    headers : list
        List of headers.
    sequences : list
        List of sequences which all have the same length.
    d_thresh : int, optional
        Maximum inclusive Hamming distance for which absorbtion can occur.
    a_thresh : numpy.float32, optional
        Minimum inclusive ratio of abundances for which absorbtion can occur.
    debug : bool, optional
        Bool used to control whether or not to return parent_children dictionary.

    Returns
    -------
    updated_headers : list
        List of headers updated with merged abundances.
    updated_sequences : list
        List of sequences corresponding to the aforementioned headers.
    parent_child : dict
        Dictionary containing which headers got absorbed into what other headers.
    """

    if a_thresh is None:
        a_thresh = 1.0

    sorted_headers, sorted_sequences = sort_data_by_abundance(headers, sequences)
    sorted_abundances = [get_abundance(u) for u in sorted_headers]

    parent_child = {}
    for i,seq1 in enumerate(sorted_sequences):
        #  The sequence has already been absorbed.
        if sorted_abundances[i] == 0:
            continue
        parent_name = sorted_headers[i].split("|")[0]
        for j,seq2 in reversed(list(enumerate(sorted_sequences))):
            if i == j:
                break
            if sorted_abundances[j] == 0:
                continue
            abundance_ratio = sorted_abundances[i] / sorted_abundances[j]
            #  The abundance ratio must be greater than or equal to a_thresh.
            if abundance_ratio < a_thresh:
                break
            #  The Hamming distance must be less than or equal to d_thresh.
            if hamming_distance(seq1, seq2) <= d_thresh:
                if parent_name not in parent_child:
                    parent_child[parent_name] = []
                #  Add merged sequences as children to parent sorted_headers[i].
                parent_child[parent_name].append(sorted_headers[j].split("|")[0])
                #  Update the abundances.
                sorted_abundances[i] += sorted_abundances[j]
                sorted_abundances[j] = 0

    #  Get resulting sequences and headers 
    updated_sequences = []
    updated_headers = []

    for index, abun in enumerate(sorted_abundances):
        #    continue
        updated_sequences.append(sorted_sequences[index])
        updated_headers.append(update_header(sorted_headers[index], sorted_abundances[index]))
    if debug:
        return updated_headers, updated_sequences, parent_child
    else:
        return updated_headers, updated_sequences

def group_data(headers: List[str], sequences: List[str]) -> dict:
    """Splits lists of headers and sequences in groupings of similar header information.

    At the very least, this function splits headers/sequenes into groupings based
    on C primer, V primer, and sequence length. At most, this function splits
    headers/sequences into groupings based on C primer, V primer, sequence length,
    time, and replicate. Error correction should be performed exclusively on sequences
    with from the same primers, times, and replicates. The length requirement comes from
    that fact that the Hamming distance is used to perform error correction.

    Parameters
    ----------
    lineage : (dict)
              nested dictionaries of lineages [V][J][L][cluster_id]

    Returns
    -------
    condensed_lineages : (dict)
                         dictionary of lineages by [(V, J, L, cluster_id)]
    """
    #  Group by sequence length, c primer, v primer,
    #  time, and replicate. What grouping is determined
    #  by how much information is in the header.
    num_header_info = len(headers[0].split("|"))
    if num_header_info >= 7:
        key_func=lambda h_in,s_in: (len(s_in), get_cprimer(h_in),
                                    get_vprimer(h_in), get_time(h_in),
                                    get_replicate(h_in))
    elif num_header_info >= 5:
        key_func=lambda h_in,s_in: (len(s_in), get_cprimer(h_in),
                                    get_vprimer(h_in), get_time(h_in))
    elif num_header_info == 4:
        key_func=lambda h_in,s_in: (len(s_in), get_cprimer(h_in),
                                    get_vprimer(h_in))
    elif num_header_info < 4:
        print("Missing information in header.\n",
              "Cannot group sequences correctly." )
        return

    grouped_data = {}
    for h,s in zip(headers, sequences):
        #  Get appropriate header information for grouping.
        key = key_func(h,s)
        if key not in grouped_data:
            grouped_data[key] = {'headers': [], 'sequences': []}
        grouped_data[key]['headers'].append(h)
        grouped_data[key]['sequences'].append(s)

    return grouped_data

def main():
    """
    usage: python error_correct.py -h"""
    import argparse

    parser = argparse.ArgumentParser(
        description='FASTA error correction')
    parser.add_argument('--fasta', type=str, help='path to FASTA')
    parser.add_argument('--outfile', type=str, help='path and name of output file')
    parser.add_argument('--delta_r', type=float, default=1,
                        help='marginal Hamming distance tolerance per decade in log'
                        'ratio abundances (default 1)')
    parser.add_argument('--delta_a', type=float, default=None,
                        help='marginal abundance tolerance of clusterable sequences per'
                        'decade in log ratio abundances (default 1)')
    parser.add_argument('--d_thresh', type=float, default=1,
                        help='Hamming distance tolerance (default 1)')
    parser.add_argument('--a_thresh', type=float, default=None,
                        help='abundance ratio tolerance (if None, not '
                             'applied)')
    parser.add_argument('--passes', type=int, default=1,
                        help='number of times to repeat greedy clustering '
                             '(default 1)')
    args = parser.parse_args()

    headers, sequences = fasta_parse(args.fasta)
    grouped_data = group_data(headers, sequences)

    if args.delta_r is not None:
        delta_r = args.delta_r
    else:
        delta_r = 1.0
    if args.delta_a is not None:
        delta_a = args.delta_a
    else:
        delta_a = 1.0
    if args.d_thresh is not None:
        d_thresh = args.d_thresh
    else:
        d_thresh = 2
    if args.a_thresh is not None:
        a_thresh = args.a_tresh
    else:
        a_thresh = None

    print('initial_unique_counts', len(headers))
    print('initial_abundance_counts', sum([get_abundance(h)
                                           for h in headers]))
    marginal_headers = []
    for this_pass in range(args.passes):
        for key in grouped_data:
            marginal_output = error_correct_marginal(grouped_data[key]['headers'],
                                                     grouped_data[key]['sequences'],
                                                     delta_r=delta_r,
                                                     delta_a=delta_a)
            if this_pass == args.passes - 1:
                marginal_headers += marginal_output[0]
            total_output = error_correct_total(marginal_output[0],
                                               marginal_output[1],
                                               d_thresh=d_thresh,
                                               a_thresh=a_thresh)

            #  Update dictionary so things don't have to be resorted
            grouped_data[key]['headers'] = total_output[0]
            grouped_data[key]['sequences'] = total_output[1]

    print('marginal_unique_counts', sum([1
                                         for h in marginal_headers
                                         if get_abundance(h) > 0]))
    print('marginal_abundance_counts', sum([get_abundance(h)
                                            for h in marginal_headers]))

    #  Combine all output and remove headers/sequences with 0 abundance.
    out_headers = []
    out_sequences = []
    for key in grouped_data:
        for i,h in enumerate(grouped_data[key]['headers']):
            if get_abundance(h) > 0:
                out_headers.append(h)
                out_sequences.append(grouped_data[key]['sequences'][i])
    print('final_unique_counts', sum([1
                                      for h in out_headers
                                      if get_abundance(h) > 0]))
    print('final_abundance_counts', sum([get_abundance(h)
                                            for h in out_headers]))

    write_to_fasta(args.outfile, out_headers, out_sequences)

if __name__ == '__main__':
    main()
