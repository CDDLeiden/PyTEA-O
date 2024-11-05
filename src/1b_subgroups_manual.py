#!/usr/bin/env python3

import os
import sys
import argparse
from os import makedirs
from Bio import SeqIO
from Bio import AlignIO

"""
Manually define subgroups

"""

def search_id(proteins_ids:dict,id: str) -> str:
    """
    Search for an ID in a dictionary of format d[protein] = [id1, id2, etc.] 
    and return the protein where the ID belongs to.

    Args:
        proteins_ids (dict): d[protein] = [id1, id2, etc.]
        id (str): Sequence ID to search for.

    Returns:
        str: The protein associated with the given ID.

    Raises:
        SystemExit: If the ID is not found in any of the protein IDs.
    """  
    for protein, ids in proteins_ids.items():
        if id in ids:
            return protein
    
    # If the ID was not found, print an error and exit
    print(f"{id} not found in any of the IDs of the different input files")
    exit()


def parse_sequence_files(fastafiles:list) -> dict:
    """
    Read multiple input sequence files and store IDs in a dictionary 
    where d[protein] = [ids].

    Args:
        sequence_files (list): List of sequence file paths in FASTA format.

    Returns:
        dict: A dictionary mapping protein names to a list of sequence IDs.
    """   
    results = {}
    for file in fastafiles:
        # Get the subgroup name from the file name
        subgroup = os.path.splitext(os.path.basename(file))[0]
        print(subgroup)

        # Initialize the protein entry in the dictionary
        results[subgroup] = []

        try:
            with open(file, 'r') as handle:
                # Parse the FASTA file and collect IDs
                results[file].extend(record.id for record in SeqIO.parse(handle, "fasta"))
        except FileNotFoundError:
            print(f"Warning: The file '{file}' was not found and will be skipped.")
        except Exception as e:
            print(f"Error processing file '{file}': {e}")

    return results


def store_subgroups(file_path:str,subgroups:dict):
    #TODO
    pass


def store_output_in_txt_file(fp, ids, y):
    with open(fp, 'w+') as file:
        for i, j in zip(ids, y):
            file.write(i + '\t' + j + '\n')

def main(args):

    # Read sequence files and store subgroups with their corresponding IDs in dictionary
    subgroups = parse_sequence_files(args.subgroupdefinitionfiles)
    
    print(subgroups)

    sys.exit()

    # Read msa containing sequencing of multiple proteins
    msa = AlignIO.read(args.msapath, "fasta")

    # For every ID in the MSA find to which protein it belongs to
    y, ids = [], []
    for i in msa:
        ids.append(i.id)
        group = search_id(subgroups, i.id)
        y.append(group)

    # Create file with subgroups
    outpath = os.path.join(args.outpath, 'subgroups.txt')

    makedirs(args.outpath,mode=0o755,exist_ok=True)

    store_output_in_txt_file(outpath, ids, y)



if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="This script categorizes sequences in the provided MSA into specified subgroups based on the provided FASTA files.")
    parser.add_argument("-m","--msapath",required=True,help="file location of Multiple Sequence Alignment (FASTA format)")
    parser.add_argument("-o","--outpath",required=True,help="output file location")
    parser.add_argument("-s","--subgroupdefinitionfiles",required=True,nargs="+",help="FASTA file containing sequences that belong to a specific subgroup. \
                                The name of the file is equal to the name of the subgroup ([subgroupname].fasta)")
    args = parser.parse_args()

    main(args)