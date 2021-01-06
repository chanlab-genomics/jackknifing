#!/usr/bin/python2
from itertools import groupby
from D2S_tools import *
import os
import math
import logging
import argparse
import sys
import time
import itertools
from functools import wraps

import multiprocessing
from concurrent import futures

import numpy as np
DESCRIPTION = '''
Calculate the D2S score between two Kmer sets.
'''

# Pass arguments.


def unpack(target_func):
    """
    A wrapper to automatically unpack arguments into a target Callable object.
    """

    @wraps(target_func)
    def wrapper(args=None, kwargs=None):

        if args is None:
            args = tuple()

        if kwargs is None:
            kwargs = dict()

        return target_func(*args, **kwargs)

    return wrapper


def main():
    # Pass command line arguments.
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter, description=DESCRIPTION)
    parser.add_argument('--kmerset1', metavar='KmerSet1.21mers.gz', type=lambda x: check_file_exists(x),
                        required=True, help='Kmers for dataset 1, can be gziped')
    parser.add_argument('--kmerset1_freq', metavar='KmerSet1.21mers.charFreq', type=lambda x: check_file_exists(x),
                        required=True, help='Character frequency for dataset 1, can be gziped')

    parser.add_argument('--kmerset2', metavar='KmerSet2.21mers.gz',
                        type=lambda x: check_file_exists(x), required=True, help='Kmer for dataset 2, can be gziped')
    parser.add_argument('--kmerset2_freq', metavar='KmerSet2.21mers.charFreq', type=lambda x: check_file_exists(x),
                        required=True, help='Character frequency for dataset 2, can be gziped')

    parser.add_argument('--D2S_out', metavar='D2S.txt', type=lambda x: write_file_check_compression(
        x), required=False, default=sys.stdout, help='Output for D2S score (default: %(default)s)')
    parser.add_argument('--threads', type=int, required=False, default=os.cpu_count(),
                        help='The number of threads used to run the jackknife algorithm. '
                        'If 0 threads are specified then it will default to os.cpu_count().')
    parser.add_argument('--debug', action='store_true', required=False, default=False,
                        help='Print DEBUG info (default: %(default)s)')
    args = parser.parse_args()

    # Set up basic debugger
    if args.debug:
        logging.basicConfig(format='#%(levelname)s :: %(asctime)s :: %(message)s',
                            stream=sys.stdout, level=logging.DEBUG)
    else:
        logging.basicConfig(format='#%(levelname)s :: %(asctime)s :: %(message)s',
                            stream=sys.stdout, level=logging.INFO)
    logger = logging.getLogger(__name__)

    logger.debug('%s', args)  # DEBUG

    d2Score_kmerset1_VS_kmerset2 = calculate_D2S(
        args.kmerset1, args.kmerset1_freq, args.kmerset2, args.kmerset2_freq, logger, args.threads)
    logger.info('kmerset1 VS. kmerset2 d2Score:%s',
                d2Score_kmerset1_VS_kmerset2)  # INFO

    d2Score_kmerset1_VS_kmerset1 = calculate_D2S(
        args.kmerset1, args.kmerset1_freq, args.kmerset1, args.kmerset1_freq, logger, args.threads)
    logger.info('kmerset1 VS. kmerset1 d2Score:%s',
                d2Score_kmerset1_VS_kmerset2)  # INFO

    d2Score_kmerset2_VS_kmerset2 = calculate_D2S(
        args.kmerset2, args.kmerset2_freq, args.kmerset2, args.kmerset2_freq, logger, args.threads)
    logger.info('kmerset2 VS. kmerset2 d2Score:%s',
                d2Score_kmerset1_VS_kmerset2)  # INFO

    D2S_distance = d2ScoreNormalization(
        d2Score_kmerset1_VS_kmerset2, d2Score_kmerset1_VS_kmerset1, d2Score_kmerset2_VS_kmerset2)
    logger.info('D2S_distance:%s', D2S_distance)  # INFO

    args.D2S_out.write(args.kmerset1 + ';' + args.kmerset2 +
                       ';' + str(D2S_distance)+'\n')


def d2ScoreNormalization(d2Score_kmerset1_VS_kmerset2, d2Score_kmerset1_VS_kmerset1, d2Score_kmerset2_VS_kmerset2):
    '''
    Function used to normalize the D2score -> creates a distance.
    '''
    # If seld distance is zero - Not sure why this would occure.
    if d2Score_kmerset1_VS_kmerset1 == 0:
        d2Score_kmerset1_VS_kmerset1 = 0.00001
    if d2Score_kmerset2_VS_kmerset2 == 0:
        d2Score_kmerset2_VS_kmerset2 = 0.00001

    value = d2Score_kmerset1_VS_kmerset2 / \
        (math.sqrt(d2Score_kmerset1_VS_kmerset1 * d2Score_kmerset2_VS_kmerset2))
    # D2S can give negative values sometimes.
    if value <= 0:
        value = 0.00001

    D2S_distance = abs(math.log(value))
    return D2S_distance


def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]


def calculate_D2S(KmerSet1_fileName, KmerSet1_freq_fileName, KmerSet2_fileName, KmerSet2_freq_fileName, logger, num_threads):

    # Open files. Best to do this here instead of as part of argparser as we need to operate on the same file
    # at the same time - will not work if we only have one file handle.
    KmerSet1_fh = read_file_check_compression(KmerSet1_fileName)
    KmerSet2_fh = read_file_check_compression(KmerSet2_fileName)
    KmerSet1_freq_fh = read_file_check_compression(KmerSet1_freq_fileName)
    KmerSet2_freq_fh = read_file_check_compression(KmerSet2_freq_fileName)

    # Initiate the Kmer_set iterator so we can check the Kmer sizes in both sets
    Kmer_iter = iterate_Kmer_sets(KmerSet1_fh, KmerSet2_fh, logger,
                                  Both_KmerSets=True, KmerSet1_Only=False, KmerSet2_Only=False)
    KmerSet1_seq, KmerSet1_count, KmerSet2_seq, KmerSet2_count = next(
        Kmer_iter)

    # Check if the kmer_seq's are the same size.
    if len(KmerSet1_seq) != len(KmerSet2_seq):
        logger.error('Kmer sizes are different between the two datasets: %s:%s\t%s:%s',
                     KmerSet1_fh.name, len(KmerSet1_seq), KmerSet2_fh.name, len(KmerSet2_seq))  # ERROR
        sys.exit(1)

    # Save Kmer size for later
    k = len(KmerSet1_seq)
    logger.info('k-mer:%s', k)  # DEBUG

    # Load the frequencies for each dataset.
    kmerset1_freq = load_Character_Frequency(KmerSet1_freq_fh, logger)
    kmerset2_freq = load_Character_Frequency(KmerSet2_freq_fh, logger)

    # Get 'NUM_SEQUENCES' and 'NUM_CHARACTERS' from file and remove from dict.
    kmerset1_NumSeqs = kmerset1_freq.pop('NUM_SEQUENCES')
    kmerset2_NumSeqs = kmerset2_freq.pop('NUM_SEQUENCES')

    kmerset1_NumChar = kmerset1_freq.pop('NUM_CHARACTERS')
    kmerset2_NumChar = kmerset2_freq.pop('NUM_CHARACTERS')

    # Calculate the number of Kmers possible from each dataset.
    # No. K-mers = (len-k)+1
    # No. K-mers (multiple seqs) = total_bases - (num_seqs * (k-1))
    kmerset1_NumKmers = kmerset1_NumChar - (kmerset1_NumSeqs * (k-1))
    kmerset2_NumKmers = kmerset2_NumChar - (kmerset2_NumSeqs * (k-1))

    logger.debug('kmerset1_NumKmers:%s\tkmerset2_NumKmers:%s',
                 kmerset1_NumKmers, kmerset2_NumKmers)  # DEBUG

    start = time.time()
    # (Re-)Initiate the Kmer iterator. We are only interested in Kmers that are shared between both sets.
    Kmer_iter = iterate_Kmer_sets(KmerSet1_fh, KmerSet2_fh, logger,
                                  Both_KmerSets=True, KmerSet1_Only=False, KmerSet2_Only=False)

    Kmer_iter = list(chunks(list(Kmer_iter), num_threads))

    thread_args = [
        (Kmer_iter_portion, kmerset1_freq, kmerset2_freq, kmerset1_NumKmers, kmerset2_NumKmers, logger) for
        Kmer_iter_portion, kmerset1_freq, kmerset2_freq, kmerset1_NumKmers, kmerset2_NumKmers, logger in
        zip(
            Kmer_iter,
            itertools.repeat(kmerset1_freq),
            itertools.repeat(kmerset2_freq),
            itertools.repeat(kmerset1_NumKmers),
            itertools.repeat(kmerset2_NumKmers),
            itertools.repeat(logger),
        )
    ]

    d2Score = 0.0

    with futures.ThreadPoolExecutor(max_workers=num_threads) as executor:

        submitted_results = executor.map(calculate_D2S_helper, thread_args)

        for result in submitted_results:
            d2Score += result

    end = time.time()
    print("Compute dist time: ", end - start)
    # Close up
    KmerSet1_fh.close()
    KmerSet2_fh.close()
    KmerSet1_freq_fh.close()
    KmerSet2_freq_fh.close()

    return d2Score


@unpack
def calculate_D2S_helper(Kmer_iter_portion, kmerset1_freq, kmerset2_freq, kmerset1_NumKmers, kmerset2_NumKmers, logger):

    d2Score = 0.0

    for KmerSet1_seq, KmerSet1_count, KmerSet2_seq, KmerSet2_count in Kmer_iter_portion:

        # Probability of k-mer occurrence in seq 1.
        PwX = calculate_PropKmerOccurrence(KmerSet1_seq, kmerset1_freq)
        # Probability of k-mer occurrence in seq 2.
        PwY = calculate_PropKmerOccurrence(KmerSet2_seq, kmerset2_freq)

        kmerScoreXBis = KmerSet1_count - (kmerset1_NumKmers*PwX)
        kmerScoreYBis = KmerSet2_count - (kmerset2_NumKmers*PwY)

        d2Score_tmp = (kmerScoreXBis*kmerScoreYBis) / \
            math.sqrt(kmerScoreXBis*kmerScoreXBis +
                      kmerScoreYBis*kmerScoreYBis)
        d2Score += (kmerScoreXBis*kmerScoreYBis) / \
            math.sqrt(kmerScoreXBis*kmerScoreXBis +
                      kmerScoreYBis*kmerScoreYBis)

        logger.debug('PwX:%s\tPwY:%s\tkmerScoreXBis:%s\tkmerScoreYBis:%s\td2Score:%s',
                     PwX, PwY, kmerScoreXBis, kmerScoreYBis, d2Score_tmp)  # DEBUG

    return d2Score


def load_Character_Frequency(freq_fh, logger):
    '''
    Loads (from file) the frequecny of each possible character (e.g. A, T, G, C)
    observed in the orginal sequence (fasta) file.

    Expected format:
            character<\t>frequency
    '''
    freq_fh.seek(0)
    logger.debug('Loading frequencies from %s', freq_fh.name)  # DEBUG
    charFreq = {}
    for char, freq in pass_column_file(freq_fh):
        charFreq[char] = float(freq)
        logger.debug('%s:%s', char, freq)  # DEBUG

    return charFreq


def calculate_PropKmerOccurrence(Kmer_seq, charFreq):
    '''
    Calculate the probability of k-mer occurrence given the Kmer_seq and the frequency of 
    characters in the full (fasta) dataset.
    '''
    freq = 1.0
    for char in list(Kmer_seq):
        freq *= charFreq[char]
    return freq


if __name__ == '__main__':
    main()
