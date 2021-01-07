#!/usr/bin/env python3

__author__ = 'Michael Ciccotosto-Camp'
__version__ = ''

import os
import csv
import pandas as pd
from glob import glob
from pprint import pprint

from typing import Callable, Dict, Iterable, List, Optional, Set, Tuple


def get_names(target_dir: str) -> list:
    """
    Gets a list of all the gene names from the a target directory which 
    holds all the comparison files.

    Parameters:
        A directory that holds all the result file named with the following
        convention:
            "id_name_1"-"id_name_2"."ext" (DON'T include paraenthesis in
            the actual filename, they are only used here for the purpose 
            of clarity and demonstration)

    Return:
        Returns a sorted list of all the unique gene names found.
    """

    # Get a list of all the files within the target directory
    target_files = glob(os.path.join(target_dir, "**"))

    # Just get the base directory of these target files
    target_files = map(os.path.basename, target_files)

    # Remove any filename extensions
    target_files = [file_.rsplit(
        '.', maxsplit=1)[0] if '.' in file_ else file_ for file_ in target_files]

    gene_names = []

    for target_file in target_files:
        gene_names.extend(target_file.split('-'))

    return sorted(list(set(gene_names)))


def get_gene_ids_from_path(result_path: str) -> List[str]:
    """
    Extarct the gene ids from a given result path.

    Return:
        Returns the extracted gene ids as a list of strings.
    """

    result_path = os.path.basename(result_path)

    result_path = result_path.rsplit('.', maxsplit=1)[
        0] if '.' in result_path else result_path

    return result_path.split('-')


def setup_df(name_list):
    """
    Creates a blank dataframe with the appropriate row and columns names.
    """

    blank_df = pd.DataFrame(index=name_list, columns=name_list, dtype=float)
    blank_df.fillna(0.0, inplace=True)

    return blank_df


def extract_result(result_path: str):
    """
    Extracts a single result from the specified result path.

    Parameters:
        result_path:
            A path to a single file containing a distance value.

    Returns:
        Returns the distance value (as a float) from the specified result path.
    """

    with open(result_path, 'r') as result_file:
        result_str = result_file.read()

    # The result should be the the very end column after performing a split
    # on the semi-colon
    value = result_str.rsplit(';', maxsplit=1)[-1]

    try:
        return float(value)
    except ValueError:
        # print("Skipping", os.path.basename(result_path))
        #######################################################################
        # NOTE: Might want to change in the future!
        #######################################################################
        return 0


def poplate_single_result(phylip_df: pd.DataFrame, result_path: str):
    """
    Populates the dataframe with a single value from the results folder.

    Parameters:
        phylip_df:
            A dataframe that will eventually contain all the distance results.

        result_path:
            A path to a single file containing a distance value.
    """

    # Get the result value (as a float) from the specified result path
    result_value = extract_result(result_path)

    # Retrieve the gene ids from the file name
    gene_id_1, gene_id_2 = get_gene_ids_from_path(result_path)

    phylip_df.loc[gene_id_1][gene_id_2] = result_value
    phylip_df.loc[gene_id_2][gene_id_1] = result_value

    return


def populate_all_results(phylip_df: pd.DataFrame, result_dir: str):
    """
    Populates the dataframe with all results from a given result directory.

    Parameters:
        phylip_df:
            A dataframe that will eventually contain all the distance results.

        result_dir:
            A directory containing all the distance results.
    """

    # Get all the result files from the the result directory.
    target_files = glob(os.path.join(result_dir, "**"))

    for target_file in target_files:
        poplate_single_result(phylip_df, target_file)

    return


def print_phylip(phylip_df: pd.DataFrame, output_path: str):
    """
    Prints the input dataframe in standard PHYLIP format to the an output file.

    Parameter:
        phylip_df:
            A dataframe to be printed to the output directory.

        output_path:
            The output file path for the PHYLIP matrix.
    """

    # Get the number of columns from this dataframe
    rows, _ = phylip_df.shape

    # Left justify all of the indexes by 10
    phylip_df.rename(index=lambda index_: index_.ljust(10), inplace=True)

    with open(output_path, 'w', newline='') as output_file:

        # Write the nukmber of rows/cols in the first line
        print('\t' + str(rows), end='\n', flush=True, file=output_file)

        # Write the remaining matrix, omit the column (header) names
        phylip_df.to_csv(output_file, sep='\t', float_format="%.8f", header=False,
                         index=True, index_label=False)

    return


def main():

    test_dir = os.path.join('/', '90days', 's4430291',
                            'Genomes_for_AFphylogeny_D2S')
    name_list = get_names(test_dir)

    blank_df = setup_df(name_list)

    test_output_file = os.path.join(
        os.getcwd(), 'reference_mat.txt')

    # test_output_file = os.path.join(
    #     '/', '90days', 's4430291', 'reference_mat.txt')

    populate_all_results(blank_df, test_dir)

    print_phylip(blank_df, test_output_file)


if __name__ == '__main__':
    main()
