#!/usr/bin/env python3

__author__ = 'Michael Ciccotosto-Camp'
__version__ = ''

import concurrent.futures
import itertools
import os
import sys
import argparse
import numpy as np
from time import perf_counter
from concurrent import futures
from glob import glob
from functools import wraps
from threading import Lock
from pprint import pprint
from random import choices, randrange
from typing import (Any, Callable, Dict, Iterable, List, Optional, Tuple, Type,
                    Union, Set)

from Bio import SeqIO, SeqRecord

"""
Example usage:
    (Windows)
    python jackknife.py --input_path D:\2020_SS\BioInfo\jackknifing\data\example.fasta --output_path D:\2020_SS\BioInfo\jackknifing\data\example_reduce_1.fasta

    python jackknife.py --input_path D:\2020_SS\BioInfo\jackknifing\example_jar_data --output_path D:\2020_SS\BioInfo\jackknifing\example_out --portion=50

    (Unix)
    python3 ./jackknife.py --input_path ./example_jar_data --output_path ./example_jar_data_out --portion=50

    python3 ./jackknife.py --input_path ./data/example.fasta --output_path ./data/example_reduce_1.fasta -v
    
    python3 ./jackknife.py --input_path ./data/S.necroappetens_CCMP2469.genome.fasta --output_path ./data/S.necroappetens_CCMP2469.genome_red.fasta -v --threads=1
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


@unpack
def add_to_rm_dict(rm_dict: Dict[str, int], rand_keys: Iterable[str],
                   dict_keys: List[str], dict_keys_set: Set[str], weights: List[int],
                   max_dictionary: Dict[str, int], rm_dict_mutex: Lock):
    """
    A helper thread function to organise how many chunks get removed from each
    chunk.

    Parameters:
        rm_dict:
            A dictionary that describes how many chunks get removed 
            from each sequence.

        rand_keys:
            Keys to which we will attempt to remove chunks from.

        dict_keys:
            A list of keys that have not had all the chunks removed.

        dict_keys_set:
            A set of keys that have not had all the chunks removed.

        weights:
            A set of keys that have not had all the chunks removed.

        max_dictionary:
            A dictionary that indicates the maximum number of chunks that
            can be removed from each sequence.

        rm_dict_mutex:
            A mutex to remove race conditions between threads.

    Returns:
        The number of chunks added to the rm_dict.
    """

    # Count the number of chunks added to the rm_dict
    total_rm = 0

    for rand_key in rand_keys:

        dict_index = 0

        if rand_key in dict_keys_set:
            dict_index = dict_keys.index(
                rand_key)
        else:
            continue

        with rm_dict_mutex:
            rm_dict[rand_key] += 1

            if rm_dict[rand_key] >= max_dictionary[rand_key]:
                dict_keys.pop(dict_index)
                weights.pop(dict_index)
                dict_keys_set.remove(rand_key)

        total_rm += 1

    return total_rm


def create_rm_dict2(total_chunks_rm: int, portion_dictionary: Dict[str, float],
                    max_dictionary: Dict[str, int], verbose: bool, threads: int = 2):
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

        verbose:
            If true, runs the function in verbose mode.

    Return:
        A randomly generated dictionary that specifies how many chunks are to
        be removed from each sequence.
    """

    # Create a list for our dictionary keys and weights
    dict_keys: List[str] = list(portion_dictionary.keys())
    # Creating a set will make it faster to search
    dict_keys_set: Set[str] = set(dict_keys)
    weights: List[int] = list(portion_dictionary.values())

    rm_dict: Dict[str, int] = dict(zip(dict_keys, itertools.cycle([0])))

    # Create a lock for modifying the rm_dict values
    rm_dict_mutex = Lock()

    for dict_key in dict_keys:

        if rm_dict[dict_key] >= max_dictionary[dict_key]:

            dict_index = dict_keys.index(dict_key)

            dict_keys.pop(dict_index)
            weights.pop(dict_index)
            dict_keys_set.remove(dict_key)

    num_args: int = total_chunks_rm
    progress: int = 0

    progress_intervals: list = list(range(0, 100, 10))

    total_chunks_rm_count = total_chunks_rm

    if verbose:
        print("Total chunks to remove: ", total_chunks_rm_count)

    while total_chunks_rm_count > 0:

        rand_keys = np.array(choices(population=dict_keys,
                                     weights=weights, k=(total_chunks_rm_count // 10) + 1), dtype=str)
        rand_keys = np.array_split(rand_keys, threads)

        thread_args = [(rm_dict, rand_keys, dict_keys, dict_keys_set, weights, max_dictionary, rm_dict_mutex) for
                       rm_dict, rand_keys, dict_keys, dict_keys_set, weights, max_dictionary, rm_dict_mutex in
                       zip(
                           itertools.cycle([rm_dict]),
                           rand_keys,
                           itertools.cycle([dict_keys]),
                           itertools.cycle([dict_keys_set]),
                           itertools.cycle([weights]),
                           itertools.cycle([max_dictionary]),
                           itertools.cycle([rm_dict_mutex])
        )]

        with futures.ThreadPoolExecutor(threads) as executor:
            submitted_results = executor.map(add_to_rm_dict, thread_args)

            for result in submitted_results:
                removed = result

                progress += removed
                total_chunks_rm_count -= removed

        if verbose:

            progress_per = progress / num_args * 100

            while len(progress_intervals) > 0 and progress_per > progress_intervals[0]:
                print(str(progress_intervals.pop(0)) +
                      '..', flush=True, end='')

    if verbose:
        print('100', flush=True)
        print()

    return rm_dict


def generate_rm_dict(total_chunks_rm: int, portion_dictionary: Dict[str, float],
                     max_dictionary: Dict[str, int], verbose: bool, threads: int = 2):
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

    rm_dict = create_rm_dict2(
        total_chunks_rm, portion_dictionary, max_dictionary, verbose, threads=threads)

    # Check that this newly generated removal dictionary is within the bounds
    # of the max dictionary. Sometimes the removal dictionary will not be
    # within the bounds of the max dictionary since it is randomly generated
    # but this shouldn't usually be the case (especially) with a low portion
    # removal since such dictionaires probabilistically unlikely.

    while any(rm_val > max_val for rm_val, max_val in zip(rm_dict.values(), max_dictionary.values())):

        if verbose:
            print("\nRestarting rm dict creation.\n")

        rm_dict = create_rm_dict2(
            total_chunks_rm, portion_dictionary, max_dictionary, verbose, threads=threads)

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

    if verbose:
        print('Reading Data...')

    # Open the fasta file as a dictionary
    fasta_dict: Dict[str, SeqRecord.SeqRecord] = SeqIO.to_dict(
        SeqIO.parse(fasta_path, "fasta"))

    # Create a portion dictionary and max dictionary from the fasta dictionary
    total_data_len, portion_dictionary, max_dictionary = compute_fasta_stats(
        fasta_dict, chunk_size)

    # Calculate the number of chunks to remove
    total_chunks_rm = int(total_data_len * portion) // chunk_size

    if verbose:
        print('Creating remove dictionary...')

    # Create a removal dictionary delete sequences from
    rm_dict = generate_rm_dict(
        total_chunks_rm, portion_dictionary, max_dictionary, verbose, threads=threads)

    if verbose:
        print('Created remove dictionary')
        print()

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
        raise ValueError("Portion values be a value between 0 and 100.")

    if args.threads < -1:
        raise ValueError(
            "The number of threads must be strictly greater than -1.")

    # Hold a list of all the file that need to be processed with there
    # respective output files
    to_complete: List[Tuple[str, str]] = []

    if os.path.isfile(args.input_path):
        to_complete.append((args.input_path, args.output_path))

    elif os.path.isdir(args.input_path):

        output_path_dir = args.output_path

        if output_path_dir is None:
            output_path_dir = args.input_path

        if not os.path.exists(output_path_dir):
            os.makedirs(output_path_dir)

        # A list of all valid fasta file extension names
        fasta_ext = ('fasta', 'fnn', 'fna', 'faa', 'fas', 'frn')

        # Get all the fasta files from the folder
        fasta_files: List[str] = itertools.chain.from_iterable(glob(os.path.join(args.input_path, "*." + ext)) for ext in fasta_ext)

        prefix: str = "_" + str(int(args.portion))

        for file_path in fasta_files:

            path_base: str = os.path.basename(file_path)
            index: int = path_base.rfind(".")
            path_base = path_base[:index] + prefix + path_base[index:]

            to_complete.append((file_path, os.path.join(
                output_path_dir, path_base)))

    for path_in, path_out in to_complete:

        portion_remover(path_in, output_path=path_out,
                        portion=args.portion / 100, chunk_size=args.chunk_size,
                        threads=args.threads, verbose=args.verbose)

    return


def main():

    parser = argparse.ArgumentParser(description="Randomly removes a portion of "
                                     "data from a fasta file.")

    parser.add_argument('--input_path', type=str, required=True,
                        help='Option 1) The path to the fasta file that will be reduced.\n\n'
                        'Option 2) A path to a folder containing several fasta files '
                        'where every fasta file will be jackknifed. '
                        'The resulting files will have will have the portion of '
                        'data that was removed printed at the end of the file name.')
    parser.add_argument('--output_path', type=str, required=False, default=None,
                        help='A file path to output the contents of the reduced flatfile. '
                        'Default output for single files is stdout. '
                        'Default output for folders is the same folder as the input.')

    parser.add_argument('-v', '--verbose', action='count', default=0,
                        help='Include to run the script in verbose mode.'
                        ' Useful for checking progress. WARN: This may cause '
                        'the program to run slower.')
    parser.add_argument('--portion', type=float, required=False, default=40,
                        help='The portion of data to be removed. The default is a 40 percent reduction.')
    parser.add_argument('--chunk_size', type=int, required=False, default=100,
                        help='The size of the chunks that get randomly removed from sequences.'
                        ' Default is a chunk size of 100.')
    parser.add_argument('--threads', type=int, required=False, default=1,
                        help='The number of threads used to run the jackknife algorithm. '
                        'If 0 threads are specified then it will default to os.cpu_count().')

    args = parser.parse_args()
    run_jackknife(args)

    exit(0)


if __name__ == '__main__':
    main()
