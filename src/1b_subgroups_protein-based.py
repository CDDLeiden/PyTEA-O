import os
import sys
import argparse
from Bio import SeqIO
from Bio import AlignIO

"""
Subgroups are calculated based on protein

"""

def search_id(proteins_ids: dict, id: str):
    """
    Search for an ID in a dictionary of format d[protein] = [id1, id2, etc.] and return the protein where the ID belongs to

    Args:
        proteins_ids (dict):d[protein] = [id1, id2, etc.]
        id (str): sequence id

    Returns:
        result (str): protein
    """      
    found = False
    for protein, ids in proteins_ids.items():
        if id in ids:
            found = True
            result = protein
    if found:
        return result
    else:
        print(f"{id} not found in any of the IDs of the different input files")
        sys.exit()


def parse_sequence_files(sequence_files: list):
    """
    Read multiple input sequence files and store IDs in dictionary d[protein] = [ids]

    Args:
        sequence_files (list): _description_

    Returns:
        _type_: _description_
    """    
    results = {}
    for seq_file in sequence_files:
        base = os.path.basename(seq_file)
        protein = base.split('.')[0]
        results[protein] = []

        with open(seq_file) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                results[protein].append(record.id)
    return results

def store_output_in_txt_file(fp, ids, y):
    with open(fp, 'w+') as file:
        for i, j in zip(ids, y):
            file.write(i + '\t' + j + '\n')

def main(args):

    # Read sequence files and store proteins with their corresponding IDs in dictionary
    results = parse_sequence_files(args.sequence_files)

    # Read msa containing sequencing of multiple proteins
    msa = AlignIO.read(args.msainpath, "fasta")

    # For every ID in the MSA find to which protein it belongs to
    y, ids = [], []
    for i in msa:
        ids.append(i.id)
        protein = search_id(results, i.id)
        y.append(protein)

    # Create file with subgroups
    outpath = os.path.join(args.subgroupoutpath, 'subgroups.txt')
    store_output_in_txt_file(outpath, ids, y)
    print("Finished")
        
def read_arguments():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("-i", "--msainpath", help="msa file (fasta format)")
    parser.add_argument("-o", "--subgroupoutpath", help="folder where to store the txt files with subgroups")
    parser.add_argument("-s", "--sequence_files", nargs="+", help="sequence files in FASTA format. Make sure that the NAME of the sequence file is equal to [protein].fasta")
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    arguments = read_arguments()
    main(arguments)            