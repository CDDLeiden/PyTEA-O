#!/usr/bin/env python3

import os
import sys
import math
import argparse
import numpy as np
import pandas as pd
from os import makedirs
from graph_results import plot_SE

multiprocessing = None
if 'linux' in sys.platform: multiprocessing = 'multiprocessing'
elif 'darwin' in sys.platform: multiprocessing = 'multiprocess'
mp = __import__(multiprocessing)

def get_reference_indexs(ref:str,msas:pd.DataFrame) -> list:

    if ref not in msas.index:
        print(f"\n\t[W] Reference {ref} is not found in the provided MSA. Skipping to prevent errors\n\n")
        return 1

    indexs = []
    for index,item in enumerate(msas.loc[ref]):
        if item == '-':
            continue
        indexs.append(index)

    return indexs

def load_msa(msa_file:str) -> pd.DataFrame:

    """
    Read in and store MSA, checking whether length is same for all sequences present
    """

    msas = {}

    ## Cheap and dirty way to check if all MSA sizes are the same
    size_range = [float('inf'),float('-inf')]
    locus = None
    with open(msa_file,'r') as IN:

        for line in IN:
            
            ## Remove trailing spaces and new line characters
            line = line.strip()

            if line == '':
                continue
        
            ## Grab Protein_ID
            if line[0] == ">":

                locus = line[1:]

                ## Notify of duplicated Protein_ID
                if locus in msas.keys():
                
                    print(f"\n\n\t[W]  Locus {locus} has been duplicated in the MSA\n")

                continue

            ## Store the MSA string
            msas[locus] = list(line)

            size_range[0] = len(line) if len(line) < size_range[0] else size_range[0]
            size_range[1] = len(line) if len(line) > size_range[1] else size_range[1]

    if size_range[0] != size_range[1]:
        print(f"\n\n\t[E]  Sequences do not have the same length!\n\t\tMin Len: {size_range[0]}\n\t\tMax Len: {size_range[1]}")
        exit()

    msas = pd.DataFrame().from_dict(msas,orient='index')

    return msas

def load_similarity_matrix(matrix:str) -> dict:

    matrix_dir = ("/".join(__file__.split("/")[0:-2]))+"/DB_files"
    labels = []
    sim_mat = {}
    file = None

    if os.path.exists(f"{matrix_dir}/{matrix.upper()}.txt"):
        file = f"{matrix_dir}/{matrix.upper()}.txt"

    elif os.path.exists(matrix):
        file = matrix
    
    if file == None:
        print(f"\n\n\t[E]  {matrix} is not a pre-existing DB, and could not find custom DB file {matrix}\n\n")
        exit()

    with open(file,'r') as IN:
        for line in IN:
            
            line = line.strip()
            
            if line == "":
                continue

            if line[0] == "" or line[0] == "#":
                continue

            if line[0] == "!":
                labels = [x for x in line.split()[1:]]
                sim_mat = {x:{y:0 for y in labels} for x in labels}
                continue

            data = line.split()
            current_label = data[0].strip()

            for index,score in enumerate(data[1:]):
                sim_mat[current_label][labels[index]] = float(score)

        return sim_mat

def process_colummns(pd_col:pd.DataFrame) -> list:

    """
    Calculates Shannon Entropy per residue position (global entropy)

    Args:
        pd_col: A single column of a pandas dataframe containing the i-th residue
                for all sequences in the MSA

    Returns: 
        A single array containing the [index number, the shannon entropy, number of nongapped residues, gapped fraction]

    Raises:
        None
    """

    # Get the residue number for the current column
    res_num = pd_col.name
    residues = pd_col.to_numpy()

    unique,counts = np.unique(residues,return_counts=True)
    num_of_seqs = sum(counts)
    residue_counts = dict(zip(unique,counts))

    # Return Shannon Entropy
    entropy = shannon_entropy(residue_counts)

    # Get number of gaps, and gap fraction
    num_of_gaps = 0 if '-' not in residue_counts.keys() else residue_counts['-']
    num_of_non_gaps = num_of_seqs - num_of_gaps
    gapped_fraction = num_of_gaps/num_of_seqs

    # # Calculate similary using blosum62 matrix
    # if self.blosum:
    #     seq = ''.join(filtered)
    #     blosum_score = self.calculate_similarity_score(seq)
    #     if np.isnan(blosum_score):
    #         gap_fraction = np.nan
    #     linetext = [col, str(blosum_score), str(n_gaps), str(self.n_seq-n_gaps), str(gap_fraction)]
    #     self.append_line_to_file(blosum_outputfile, linetext)

    return [res_num,entropy,num_of_non_gaps,gapped_fraction,residue_counts]

def max_shannon_entropy(num_of_seqs:int) -> float:

    res_prob = (1/num_of_seqs)

    max_sh_entropy = (num_of_seqs * (res_prob * (math.log(res_prob, 2))))

    return max_sh_entropy

def shannon_entropy(res_count:dict) -> float:

    """
    Calculate Shannon's Entropy per column of an alignment (H=-\sum_{i=1}^{M} P_i\,log_2\,P_i)
    """

    num_of_seqs = sum([x for x in res_count.values()])

    entropy = 0
    
    # Number of residues in column
    for residue in res_count.keys():
        
        coeff = 1
        num_of_res = res_count[residue]

        if residue == '-':
            coeff = res_count[residue]
            num_of_res = 1

        prob = num_of_res / num_of_seqs

        entropy += (coeff * (prob * (math.log(prob, 2))))

    return entropy

def run(args:dict) -> None:

    """
    Calculate Shannon Entropy for a Multiple Sequence Alignment (MSA)

    """

    path_msa = args.msapath
    path_subfamilies = args.subfamilypath
    outpath = args.outpath
    reference_ids = args.reference_ids
    threads = args.threads
    sim_mat = args.similarity_matrix

    # Default settings
    threshold_sequences_per_subfamily = 1  # Minimalsequences per subfamily
    shannon_entropy_file = "shannon_entropy.txt"
    logo_file = "consensus_logo.txt"
    shannon_entropy_file_filtered = "shannon_entropy_filtered.txt"
    shannon_entropy_file_header = ["alignment_pos", "entropy", "n_gaps", "n_nogap", "gap_fraction"]
        
    # For calculating similarity between subgroups, based on blosum62 matrix
    similarity_file = 'similarity_blosum62.txt'
    similarity_file_filtered = "similarity_blosum62_filtered.txt"
    similarity_file_header = ["alignment_pos", "blosum62_score", "n_gaps", "n_nogap", "gap_fraction"]

    makedirs(outpath,mode=0o755,exist_ok=True)

    # # Check if there already are entropy files in the output dir > otherwise ask to remove the present files or choose another directory
    # if os.path.isfile(f"{outpath}/{shannon_entropy_file_filtered}"):
    #     sys.exit(f"Entropy file already exists in {outpath}, please choose another output directory or remove files.")

    # Load and make sure all sequences in the MSA have the same length before starting calculations
    msas = load_msa(path_msa)

    n_seq, n_res = msas.shape
    max_sh_entropy = max_shannon_entropy(n_seq)
    col_data = [msas[x] for x in range(n_res)]

    with mp.Pool(threads) as executor:

        results = executor.map(process_colummns,col_data)

    msa_buffer = len(str(n_res))
    with open(f"{outpath}/{shannon_entropy_file}",'w') as OUT:
        OUT.write(f"## MSA_position\tShannon_Entropy\tFract_Shannon_Entropy\tSequences_at_pos\tFraction_of_gaps\n")
        for res_num,sh_entropy,non_gapped_res,gapped_fraction,_ in sorted(results,key=lambda x: x[0]):
            OUT.write(f"{res_num: >{msa_buffer}d}\t{sh_entropy: >.2f}\t{sh_entropy/max_sh_entropy: >3.2f}\t{non_gapped_res: >{msa_buffer}}\t{gapped_fraction:0<3.2f}\n")

    with open(f"{outpath}/{logo_file}",'w') as OUT:
        OUT.write(f"## MSA_position\tRes:Count;\n")
        for res_num,_,_,_,residue_counts in (sorted(results,key=lambda x: x[0])):
            residues = sorted(residue_counts.keys(),key=lambda x: residue_counts[x],reverse=True)
            OUT.write(f"{res_num: >{msa_buffer}d}\t{residues[0]}:{residue_counts[residues[0]]}")
            for residue in residues[1:]:
                OUT.write(f";{residue}:{residue_counts[residue]}")
            OUT.write("\n")

    for reference in reference_ids:
        indexs = get_reference_indexs(reference,msas)
        res_buffer = len(str(len(indexs)))
        with open(f"{outpath}/{reference}.shannon_entropy.txt",'w') as OUT:
            OUT.write(f"## Sequence_Pos\tMSA_Pos\tResidue\tShannon_Entropy\tFract_Shannon_Entropy\tSequences_at_pos\tFraction_of_gaps\n")
            for res_index,msa_index in enumerate(indexs):
                _,sh_entropy,non_gapped_res,gapped_fraction,_ = results[msa_index]
                OUT.write(f"{res_index: <{res_buffer}d}\t{msa_index: >{msa_buffer}d}\t{msas.loc[reference][msa_index]}\t{sh_entropy: >.2f}\t{sh_entropy/max_sh_entropy: >3.2f}\t{non_gapped_res: >{msa_buffer}}\t{gapped_fraction:0<3.2f}\n")

    # load_similarity_matrix(matrix=sim_mat)

# def calculate_similarity_score(seq, method='BLOSUM62'):
#     """
#     Using BLOSUM62 matrix. Excluded diagonal and divide by length sequence
#     """
#     aligner = Align.PairwiseAligner()
#     aligner.substitution_matrix = substitution_matrices.load(method)

#     amino_acids = SE.natural_amino_acids()

#     # Remove gaps from sequence
#     # Remove non AA characters from sequence
#     full_length = len(seq)
#     seq = seq.replace('-', '')
#     for char in seq:
#         if char not in amino_acids:
#             seq = seq.replace(char, '')
            
#     num_gaps = full_length - len(seq)
    
#     # Set gap penalty
#     gap_penalty = 0.5 / full_length  # Penalty for each gap. 0.5 is max similarity score

#     # Make matrix
#     M = np.empty((len(seq), len(seq)), dtype='int')

#     # Get upper/lower triangle (including diagonal) and calculate score
#     for i in range(len(seq)):
#         for num, j in enumerate(seq):
#             score = aligner.substitution_matrix[str(seq[i]), str(seq[num])]
#             M[i, num] = float(score)
    
#     # Calculate score by summing values
#     sum_tri = np.sum(M[np.triu_indices(len(seq), k=1)])  # Exclude diagonal
#     sum_diagonal = np.trace(M)  
#     if len(seq) == 1:  # Cannot (and should not) be calculated (Runtimerror)
#         score = np.nan
#     else:
#         score = sum_tri / sum_diagonal  # Divide by diagonal for normalization between residues
#         score = score - (gap_penalty * num_gaps)  # Subtract gap penalty
#         score = round(score / full_length, 2)
    
#     # If score is negative, set to 0
#     if score < 0:
#         score = 0

#     return score


def main(args):

    run(args)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='')
    parser.add_argument("-m", "--msapath", required=True, help="file location of Multiple Sequence Alignment (FASTA format)")
    parser.add_argument("-s", "--subfamilypath", required=True, help="file location of subfamily assignment (TXT file)")
    parser.add_argument("-o", "--outpath", required=True, help="output file location")
    parser.add_argument("-r", "--reference_ids", required=False, help="reference ID(s).", nargs="+")
    # parser.add_argument("-i", "--seqidentity", help="Calculate pairwise sequence identity matrix", required=False, default=False)
    parser.add_argument("-x", "--similarity_matrix", help="Calculate scores according to provided matrix", required=False, type=str,default=None)
    # parser.add_argument("-j", "--json", type=str, help="store command line argument as json. Provide folder where to store", required=False)
    parser.add_argument("-t","--threads",type=int,required=False,default=1)
    args = parser.parse_args()

    main(args)