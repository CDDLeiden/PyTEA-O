#!/usr/bin/env python3

import os
import json
import pickle
import subprocess
import numpy as np
from Bio import AlignIO
from Bio import pairwise2 as pw2

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

def load_pickle(in_fp):
    obj = pickle.load(open(in_fp, "rb" ))
    return obj

def rename(fp, extension):
    pre, _ = os.path.splitext(fp)
    renamed = pre + extension
    return renamed

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

def sequence_identity(seq1: str, seq2: str) -> float:
    global_align = pw2.align.globalxx(seq1, seq2)
    matches = global_align[0][2]
    seq_length = len(global_align[0][0])
    percent_match = (matches / seq_length) * 100
    #print(percent_match)
    return percent_match