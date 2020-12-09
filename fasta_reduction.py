#!/usr/bin/env python3

__author__ = 'Michael Ciccotosto-Camp'
__version__ = ''

import concurrent.futures
import itertools
import os
import sys
from array import array
from concurrent import futures
from functools import wraps
from threading import Lock
from pprint import pprint
from random import choices, randrange
from typing import (Any, Callable, Dict, Iterable, List, Optional, Tuple, Type,
                    Union)

from Bio import SeqIO, SeqRecord


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
def remove_chunks(fasta_dict, seq_id, chunk_size, num_chunks_rm, mutex):
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
                    portion: float = 0.4, chunk_size: int = 100, threads: int = 1):
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
    """

    # TODO: Maybe move this to when the cmd line arguments are being retrieved.
    if not 0 < portion < 1:
        raise ValueError("Portion values be a value between 0 and 1.")

    if threads < -1:
        raise ValueError("Portion values be a value between 0 and 1.")

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

    for thread_arg in thread_args:
        remove_chunks(thread_arg)

    with futures.ThreadPoolExecutor(threads) as executor:
        submitted_results = executor.map(remove_chunks, thread_args)

        for result in submitted_results:
            _ = result

    record_list = []

    for rec_id, rec_val in zip(fasta_dict.keys(), fasta_dict.values()):
        rec_val.id = rec_id
        record_list.append(rec_val)

    SeqIO.write(record_list, output_path, 'fasta')

    return


def main():

    example_fasta_path = os.path.join(
        os.getcwd(), 'data', 'example.fasta')

    example_out_path = os.path.join(
        os.getcwd(), 'data', 'example_red_1.fasta')

    portion_remover(example_fasta_path,
                    output_path=example_out_path, threads=2)


if __name__ == '__main__':
    main()
