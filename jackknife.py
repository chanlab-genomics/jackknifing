#!/usr/bin/env python3

__author__ = 'Michael Ciccotosto-Camp'
__version__ = ''

import concurrent.futures
import itertools
import os
import sys
import argparse
from array import array
from concurrent import futures
from functools import wraps
from threading import Lock
from pprint import pprint
from random import choices, randrange
from typing import (Any, Callable, Dict, Iterable, List, Optional, Tuple, Type,
                    Union)

from Bio import SeqIO, SeqRecord

"""
Example usage:
    (Windows)
    python jackknife.py --input_path D:\2020_SS\BioInfo\jackknifing\data\example.fasta --output_path D:\2020_SS\BioInfo\jackknifing\data\example_reduce_1.fasta

    (Unix)
    python3 python ./jackknife.py --input_path ./data/example.fasta --outputpath ./data/example_reduce_1.fasta
"""


def unpack(orig_func):
    """
    A wrapper to automatically unpack arguments into a target function.
    """

    @wraps(orig_func)
    def wrapper(args=None, kwargs=None):

        if args is None:
            args = tuple()

        if kwargs is None:
            kwargs = dict()

        return orig_func(*args, **kwargs)

    return wrapper


def compute_fasta_stats(fasta_dict: Dict[str, SeqRecord.SeqRecord], chunk_size: int) -> Tuple[int, Dict[str, float], Dict[str, int]]:
    """
    Creates a portion dictionary and a max dictionary from a fasta dictionary.

    NOTE:
        The the statistics are computed in the same function for computational
        efficiency.

    Parameters:
        fasta_dict:
            The fasta file formatted as a dictionary (by BioPython).

        chunk_size:
            The chunk size of the data to remove.

    Returns:
        A tuple containing the portion dictionary and max dictionary. ie:
            portion_dictionary, max_dictionary
    """

    portion_dictionary: Dict[str, float] = {}
    max_dictionary: Dict[str, int] = {}

    total_data_len = 0

    for seq_id, record in zip(fasta_dict.keys(), fasta_dict.values()):

        # Get the length of the current sequence
        seq_len = len(record.seq)

        total_data_len += seq_len

        portion_dictionary[seq_id] = seq_len
        max_dictionary[seq_id] = (seq_len - 1) // chunk_size

    # Now compute the portion of each sequence (instead of just the total
    # sequence length)
    for seq_id in fasta_dict.keys():

        portion_dictionary[seq_id] = portion_dictionary[seq_id] / \
            total_data_len

    return total_data_len, portion_dictionary, max_dictionary


def create_rm_dict(total_chunks_rm: int, portion_dictionary: Dict[str, float]):
    """
    Creates a dictionary that specifies how many chunks are to be removed from
    each sequence. The dictionary returned by this function will not
    nesseccarily be within the constraints of the max_dictionary.

    Parameters:
        total_chunks_rm:
            The total number of chunks that need to be removed from the fasta
            file.

        portion_dictionary:
            A dictionary indicating what portion of the data each sequence
            takes up.

        max_dictionary:
            A dictionary that indicates the maximum number of chunks that
            can be removed from each sequence.

    Return:
        A randomly generated dictionary that specifies how many chunks are to
        be removed from each sequence.
    """

    # Create a list for our dictionary keys and weights
    dict_keys = list(portion_dictionary.keys())
    weights = list(portion_dictionary.values())

    rm_dict = dict(zip(dict_keys, itertools.cycle([0])))

    # Randomly pick dictionary keys to remove chunks from using the above weights.
    # NOTE: The 'choices' function picks elements from the population with
    # replacement
    rand_keys = choices(population=dict_keys,
                        weights=weights, k=total_chunks_rm)

    for rand_key in rand_keys:
        rm_dict[rand_key] += 1

    return rm_dict


def generate_rm_dict(total_chunks_rm: int, portion_dictionary: Dict[str, float],
                     max_dictionary: Dict[str, int]):
    """
    Creates a dictionary that specifies how many chunks are to be removed from
    each sequence. The dictionary returned by this function will be within the
    constraints of the max_dictionary.

    Parameters:
        total_chunks_rm:
            The total number of chunks that need to be removed from the fasta
            file.

        portion_dictionary:
            A dictionary indicating what portion of the data each sequence
            takes up.

        max_dictionary:
            A dictionary that indicates the maximum number of chunks that
            can be removed from each sequence.

    Return:
        A randomly generated dictionary that specifies how many chunks are to
        be removed from each sequence.
    """

    rm_dict = create_rm_dict(
        total_chunks_rm, portion_dictionary)

    # Check that this newly generated removal dictionary is within the bounds
    # of the max dictionary. Sometimes the removal dictionary will not be
    # within the bounds of the max dictionary since it is randomly generated
    # but this shouldn't usually be the case (especially) with a low portion
    # removal since such dictionaires probabilistically unlikely.

    while any(rm_val > max_val for rm_val, max_val in zip(rm_dict.values(), max_dictionary.values())):
        rm_dict = create_rm_dict(
            total_chunks_rm, portion_dictionary)

    return rm_dict


@unpack
def remove_chunks(fasta_dict: dict, seq_id: str, chunk_size: int, num_chunks_rm: int, mutex: Lock):
    """
    Removes a prescribed number of chunks from a sequence.

    Parameter:
        fasta_dict:
            The fasta file formatted as a dictionary (by BioPython).

        seq_id:
            The ID of the sequence we are removing from.

        chunk_size:
            The size of the chunk to remove.

        num_chunks_rm:
            The number of chunks to remove.

        mutex:
            A threading mutex to safely access the fasta dictionary (which
            will be shared between multiple threads).
    """

    new_seq = ""

    with mutex:
        # Convert the BioPython sequence into an array
        seq_array = fasta_dict[seq_id].seq._data

    if chunk_size * num_chunks_rm == len(seq_array):

        # The very unlikely event where the entire sequence should be removed.
        new_seq = ""

    else:

        # Get a condensed list of indexes to start removal
        rm_start_indexes = list(range(len(seq_array) // chunk_size))

        # Randomly choose a subset of removal indices
        rm_start_indexes = sorted(choices(rm_start_indexes, k=num_chunks_rm))

        rm_indexes = []

        for index in rm_start_indexes:
            rm_indexes.append(index * chunk_size)
            rm_indexes.append((index + 1) * chunk_size)

        rm_indexes = [0] + rm_indexes + [len(seq_array)]

        rm_indexes = zip(rm_indexes[::2], rm_indexes[1::2])

        for rm_start, rm_end in rm_indexes:

            new_seq += seq_array[rm_start:rm_end]

    with mutex:
        fasta_dict[seq_id].seq._data = new_seq

    return


@unpack
def remove_chunks2(fasta_dict, seq_id, chunk_size, num_chunks_rm, mutex):
    """
    Removes a prescribed number of chunks from a sequence.

    Parameter:
        fasta_dict:
            The fasta file formatted as a dictionary (by BioPython).
    """

    with mutex:
        # Convert the BioPython sequence into an array
        seq_array = fasta_dict[seq_id].seq._data

    if chunk_size * num_chunks_rm == len(seq_array):

        # The very unlikely event where the entire sequence should be removed.
        seq_array = ""

    else:

        for _ in range(num_chunks_rm):

            # Pick a starting value to remove the chunk
            index_start = randrange(0, len(seq_array) - chunk_size)

            seq_array = seq_array[:index_start] + \
                seq_array[index_start + chunk_size:]

    with mutex:
        fasta_dict[seq_id].seq._data = seq_array

    return


def portion_remover(fasta_path: str, output_path: str = None,
                    portion: float = 0.4, chunk_size: int = 100, threads: int = 1, verbose: bool = True):
    """
    Randomly removes a certain portion of data from a fasta file and saves the
    result in an output path.

    Parameters:
        fasta_path:
            A path to the fasta file to remove data.

        output_path:
            The output path to save the resulting data. If no path in
            specified than the output will be directed to stdout.

        portion:
            The portion of data to remove. The default is a 40% reduction
            (meaning 60% of the data will remain).

        chunk_size:
            The chunk size of the data to remove.

        threads:
            The number of threads used to remove the data.

        verbose:
            If true, runs the function in verbose mode.
    """

    if threads == 0:
        threads = os.cpu_count()

    if output_path is None:
        output_path = sys.stdout

    if verbose:
        print()
        print('Running on ' + os.path.basename(fasta_path) + ' with:')
        print('\t' + 'output_path=' + output_path)
        print('\t' + 'portion=' + str(portion * 100))
        print('\t' + 'chunk_size=' + str(chunk_size))
        print('\t' + 'threads=' + str(threads))
        print()

    # Create a lock for multi-threading later on
    mutex = Lock()

    # Open the fasta file as a dictionary
    fasta_dict: Dict[str, SeqRecord.SeqRecord] = SeqIO.to_dict(
        SeqIO.parse(fasta_path, "fasta"))

    # Create a portion dictionary and max dictionary from the fasta dictionary
    total_data_len, portion_dictionary, max_dictionary = compute_fasta_stats(
        fasta_dict, chunk_size)

    # Calculate the number of chunks to remove
    total_chunks_rm = int(total_data_len * portion) // chunk_size

    # Create a removal dictionary delete sequences from
    rm_dict = generate_rm_dict(
        total_chunks_rm, portion_dictionary, max_dictionary)

    thread_args = [(fasta_dict, seq_id, chunk_size, num_chunks, mutex) for
                   fasta_dict, seq_id, chunk_size, num_chunks, mutex in
                   zip(itertools.cycle([fasta_dict]), rm_dict.keys(), itertools.cycle([chunk_size]), rm_dict.values(), itertools.cycle([mutex]))]

    num_args: int = len(thread_args)
    progress: int = 0

    progress_intervals: list = list(range(0, 100, 10))

    if verbose:
        print('Progress: ', flush=True, end='')
        print()

    with futures.ThreadPoolExecutor(threads) as executor:
        submitted_results = executor.map(remove_chunks, thread_args)

        for result in submitted_results:
            _ = result

            if verbose:
                progress += 1

                progress_per = progress / num_args * 100

                while len(progress_intervals) > 0 and progress_per > progress_intervals[0]:
                    print(str(progress_intervals.pop(0)) +
                          '..', flush=True, end='')
    if verbose:
        print('100', flush=True)
        print()

    record_list = []

    for rec_id, rec_val in zip(fasta_dict.keys(), fasta_dict.values()):
        rec_val.id = rec_id
        record_list.append(rec_val)

    SeqIO.write(record_list, output_path, 'fasta')

    return


def run_jackknife(args):

    if not 0 < args.portion < 100:
        raise ValueError("Portion values be a value between 0 and 1.")

    if args.threads < -1:
        raise ValueError(
            "The number of threads must be strictly greater than -1.")

    portion_remover(args.input_path, output_path=args.output_path,
                    portion=args.portion / 100, chunk_size=args.chunk_size,
                    threads=args.threads, verbose=args.verbose)

    return


def main():

    parser = argparse.ArgumentParser(description="Randomly removes a portion of "
                                     "data from a fasta file.")

    parser.add_argument('--input_path', type=str, required=True,
                        help='The path to the fasta file that will be reduced.')
    parser.add_argument('--output_path', type=str, required=False, default=None,
                        help='A file path to output the contents of the reduced flatfile. '
                        'Default output file is stdout.')

    parser.add_argument('-v', '--verbose', action='count', default=0,
                        help='Include to run the script in verbose mode.'
                        ' Useful for checking progress.')
    parser.add_argument('--portion', type=float, required=False, default=40,
                        help='The portion of data to removed. The default is 40 (40\% removal).')
    parser.add_argument('--chunk_size', type=int, required=False, default=100,
                        help='The size of the chunks that get randomly removed from sequences.'
                        ' Default is 100.')
    parser.add_argument('--threads', type=int, required=False, default=1,
                        help='The number of threads used to run the jackknife algorithm.'
                        'If 0 threads are specified then it will default to os.cpu_count().')

    args = parser.parse_args()
    run_jackknife(args)

    exit(0)


if __name__ == '__main__':
    main()
