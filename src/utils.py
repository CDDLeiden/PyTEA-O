#!/usr/bin/env python3

import os
import json
import pickle
import subprocess
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio import Align
from Bio import AlignIO
from Bio import pairwise2 as pw2
from Bio.Align import substitution_matrices

def dict_to_json(obj, out_fp):
    with open(out_fp, 'w') as fp:
        json.dump(obj, fp)

def save_commandline_args(script, args, folder):
        script_base = script.split('.')[0]
        fp = os.path.join(folder, script_base + '.json')
        args['commit'] = fetch_git_commit()  # Add git commit to dict
        dict_to_json(args, fp)

def fetch_git_commit():
    hash = subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD']).decode('ascii').strip()
    return hash

def save_pickle(obj, out_fp):
    with open(out_fp, 'wb') as f:
        pickle.dump(obj, f)

# def load_pickle(in_fp):
#     file = open(in_fp,'rb')
#     obj = pickle.load(file)
#     return obj

def load_pickle(in_fp):
    obj = pickle.load(open(in_fp, "rb" ))
    return obj

def read_sequence(fp, extension="fasta"):
    with open(fp) as handle:
        for record in SeqIO.parse(handle, extension):
            pass
    return record.id, record.seq

def sequence_records(msafile, extension="fasta"):
    return list(SeqIO.parse(msafile, extension))

def read_alignment(fp, extension="fasta"):
    return AlignIO.read(fp, extension)

def rename(fp, extension):
    pre, _ = os.path.splitext(fp)
    renamed = pre + extension
    return renamed

def extract_from_tsf(fp, col_num, type=float, begin=0):
    res = []
    with open(fp, 'r') as f:
        content = f.readlines()[begin:]
        for line in content:
            line = line.split()
            res.append(type(line[col_num]))
    return res

def lists_to_tsf(*lists, fp, fname):
    with open(os.path.join(fp, fname), 'w+') as out:
        for items in zip(*lists):
            out.write('\t'.join(items) + '\n')

def print_pymol_selection(df, column: str) -> None:
    """Print residues from pandas dataframe '+' separated so selection can directly be copied into PyMol

    Args:
        df (pd.Dataframe): dataframe with residue ID column in it
        column (str): column containing the residue IDs
    """
    fullresidues = df[column].tolist()
    filtered = [i for i in fullresidues if i != '-']
    resids = [int(i[1:]) for i in filtered]
    ic_pymol = ('+'.join([str(i) for i in resids]))
    print(ic_pymol)
    return None

def sequence_identity(seq1: str, seq2: str) -> float:
    """Calculate sequence identity of two sequences

    Args:
        seq1 (str): one of the two sequences to calculate the identity for
        seq2 (str): one of the two sequences to calculate the identity for

    Returns:
        float: percent identity between seq1 and seq2
    """
    global_align = pw2.align.globalxx(seq1, seq2)
    matches = global_align[0][2]
    seq_length = len(global_align[0][0])
    percent_match = (matches / seq_length) * 100
    return percent_match

def pairwise_seqidentity_matrix(msa_fp: str) -> np.array:
    """
    Calculate the pairwise sequence identity between all sequences from a MSA (.fasta format) and return as as array.
    Gaps are removed from the MSA and the identities are stored in the bottom half of a np.array

    The identitiy array is stored in the same location as the MSA file with the extension '_seqidentity.pkl'
    The labels of the identity array are stored in the same location as the MSA file with the extension '_seqidentity_labels.pkl'

    Args:
        msa_fp (str): filepath of the MSA

    Returns:
        np.array: _description_
    """
    identity_matrix_file_extension = '_seqidentity.pkl'
    identity_matrix_labels_file_extension = '_seqidentity_labels.pkl'

    alignment = AlignIO.read(msa_fp, "fasta")  # Read alignment
    labels = [rec.id for rec in alignment]

    M = np.empty((len(alignment), len(alignment)))  # Create output matrix

    for i in range(len(alignment)):
        print(i)
        for num, j in enumerate(alignment):
            if i >= num:  # Only calculate lower triangle, including diagonal, to save time
                identity = sequence_identity(alignment[i].seq.replace('-', '') , j.seq.replace('-', ''))  # Remove gaps from sequence
                identity = round(identity, 2)
                M[i, num] = float(identity)
            else:
                M[i, num] = 0
    
    save_pickle(M, rename(msa_fp, identity_matrix_file_extension))  # Save matrix  and labels as pkl file
    save_pickle(labels, rename(msa_fp, identity_matrix_labels_file_extension))
    return M

def blosum62_scores_matrix(msa_fp: str) -> np.array:
    """
    Calculate the BLOSUM62 scores between all sequences of an MSA and return the matrix
    """
    matrix_file_extension = '_seqsimilarity.pkl'
    matrix_labels_file_extension = '_seqsimilarity_labels.pkl'

    alignment = AlignIO.read(msa_fp, "fasta")  # Read alignment
    labels = [rec.id for rec in alignment]

    M = np.empty((len(alignment), len(alignment)))  # Create output matrix

    for i in range(len(alignment)):
        for num, j in enumerate(alignment):
            score = score_pairwise(alignment[i].seq.replace('-', '') , j.seq.replace('-', ''))  # Remove gaps from sequence
            score = round(score, 2)
            M[i, num] = float(score)
    
    save_pickle(M, rename(msa_fp, matrix_file_extension))  # Save matrix  and labels as pkl file
    save_pickle(labels, rename(msa_fp, matrix_labels_file_extension))

def score_pairwise(seq1: str, seq2: str) -> int:
    """
    Calculate the similarity score using the BLOSUM62 matrix

    Args:
        seq1 (str): _description_
        seq2 (str): _description_

    Returns:
        int: _description_
    """
    # "BLOSUM: A block-oriented substitution matrix for proteins" suggests 11 for gap opening and extension
    aligner = Align.PairwiseAligner()
    # aligner.open_gap_score = 11
    # aligner.extend_gap_score = 11
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    # print("open gap score: ", aligner.open_gap_score)
    # print("extend gap score ", aligner.extend_gap_score)
    alignments = aligner.align(seq1, seq2)
    optimal = next(alignments)
    score = optimal.score
    # print(score)
    return score

def molecular_weight_aas():
    # source https://www.sigmaaldrich.com/NL/en/technical-documents/technical-article/protein-biology/protein-structural-analysis/amino-acid-reference-chart
    d = {"A": 89.10, "R": 174.2, "N": 132.12, "D": 133.11, "C": 121.16, "E": 147.13, "Q":146.15, "G":75.07, "H": 155.16,
         "I": 131.18, "L": 131.18, "K": 146.19, "M": 149.21, "F": 165.19, "P": 115.13, "S": 105.09, "T": 119.12, "W": 204.23,
         "Y": 181.19, "V": 117.15}
    return d

def residue_groups():
    d = {"positive": ["R", "H", "K"],
         "negative": ["D", "E"],
         "polar_uncharged": ["S", "T", "N", "Q"],
         "special": ["C", "G", "P"],
         "hydrophobic_noring": ["A", "V", "I", "L", "M"],
         "hydrophobic_ring": ["F", "Y", "W"]}
    return d

def sequence_identity(seq1: str, seq2: str) -> float:
    global_align = pw2.align.globalxx(seq1, seq2)
    matches = global_align[0][2]
    seq_length = len(global_align[0][0])
    percent_match = (matches / seq_length) * 100
    #print(percent_match)
    return percent_match


def convert_aa_format():
    d = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
         'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
         'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
         'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'}
    return d

# def pairwise_seqidentity_matrix(msa_fp: str) -> np.array:
#     """Calculate pairwise seqidentity for all sequences in an MSA and return as as array."""
#     identity_matrix_file_extension = '_seqidentity.pkl'
#     identity_matrix_labels_file_extension = '_seqidentity_labels.pkl'

#     alignment = AlignIO.read(msa_fp, "fasta")  # Read alignment
#     labels = [rec.id for rec in alignment]

#     M = np.empty((len(alignment), len(alignment)))  # Create output matrix

#     for i in range(len(alignment)):
#         for num, j in enumerate(alignment):
#             if i >= num:  # Only calculate lower triangle, including diagonal, to save time
#                 identity = sequence_identity(alignment[i].seq.replace('-', '') , j.seq.replace('-', ''))  # Remove gaps from sequence
#                 identity = round(identity, 2)
#                 M[i, num] = float(identity)
#             else:
#                 M[i, num] = 0
    
#     outmatrix 
#     save_pickle(M, rename(msa_fp, identity_matrix_file_extension))  # Save matrix  and labels as pkl file
#     save_pickle(labels, rename(msa_fp, identity_matrix_labels_file_extension))

# def mtb_lineage():
#     lineage = {
#         "genus":"Mycobacterium",
#         "phylum":"Actinobacteria",
#         "class":"Actinomycetia",
#         "order":"Corynebacteriales",
#         "family":"Mycobacteriaceae",
#         "species":"Mycobacterium"
#     }
#     return lineage