#!/usr/bin/env python3

import os
import sys
import argparse
from itertools import chain
from os import makedirs
from Bio import SeqIO
from Bio import AlignIO

"""
Subgroup definition for TEA calculations

"""


def parse_subgroups(files:list):
    """
    Read input files and store IDs in dictionary d[subgoup] = [ids]
    """    
    results = {}
    for file in files:
        base = os.path.basename(file)
        sub = base.split('.')[0]
        results[sub] = []

        with open(file) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                results[sub].append(record.id)
    return results

def store_subgroups(subgroups:dict, outfp:str):
    all_ids = list(chain.from_iterable(subgroups.values()))
    converted_keys = {key: hex(index + 1) for index, key in enumerate(subgroups.keys())}
    value_to_key_mapping = {value: converted_keys[key] for key, values in subgroups.items() for value in values}
    keys_line = ";".join(value_to_key_mapping[value] for value in all_ids) + "\n"
    with open(outfp, 'w') as file:
        first_line = "## " + ";".join(all_ids) + "\n"
        file.write(first_line)
        file.write(keys_line)


def check_input(msa, subgroups):
    """
    Find any un-assigned IDs or IDs that are assigned to multiple subgroups
    """
    msa = AlignIO.read(msa, "fasta")
    print("Length of MSA: ", len(msa))
    print("Number of subgroups: ", len(subgroups))
    assigned = []
    for i in msa:
        for key, value in subgroups.items():
            if i.id in value:
                assigned.append(i.id)
    unassigned = [i for i in msa if i.id not in assigned]
    print("Number of unassigned IDs: ", len(unassigned))
    print("Number of assigned IDs: ", len(assigned))
    if len(unassigned) > 0:
        print("Error: The following IDs are not assigned to any subgroup")
        for i in unassigned:
            print(i.id)
    if len(assigned) != len(msa):
        print("Error: Not all IDs are assigned to a subgroup")
        sys.exit(1)


def main(args):

    # Read subgroups and store IDs in dictionary
    results = parse_subgroups(args.subgroups)

    # Check if there are any un-assigned IDs
    check_input(args.msa, results)

    # Create file with subgroups
    out = os.path.splitext(args.msa)[0] + "_TEA_subgroups.txt"
    store_subgroups(results, out)
    
    
def read_arguments():
    parser = argparse.ArgumentParser(description="""Custum TEA subgroup definition. Input is a Multiple Sequence Alignment (MSA) and multiple sequence files in .fasta format.
                                                 Each sequence file contains the sequences that should be grouped together. The output is a text file containing the subfamily grouping. 
                                                 The output file is named as <msa_filename>_TEA_subgroups.txt""")
    parser.add_argument("-m", "--msa", required=True, help="msa file in .fasta format")
    parser.add_argument("-s", "--subgroups", nargs="+", help="sequence files in .fasta format. Each file contains the sequences that should be grouped together")
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    arguments = read_arguments()
    main(arguments)            