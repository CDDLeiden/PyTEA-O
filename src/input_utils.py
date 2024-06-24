import os
import sys
import time
import argparse
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio import Entrez
from multiprocessing import Pool, Manager
from utils import load_pickle, save_pickle, read_sequence, sequence_records, lists_to_tsf

def run_entrez_api(code, email="Your.Name.Here@example.org"):
    """
    Run Entrez API to search for taxonomy of record

    Args:
        code (str): sequence_id
        email (str, optional): email to provide to Entrez. Defaults to "Your.Name.Here@example.org".

    Returns:
        (int): taxonomy
    """    
    try:
        Entrez.email = email
        handle = Entrez.esummary(db="protein", id=code)
        time.sleep(0.33)
        taxonList = Entrez.read(handle)
        tax = int(taxonList[0]["TaxId"])
        handle.close()
        if tax == 1:
            print(code, tax)
    except RuntimeError:
        tax = None
    return tax

def seqid_from_fullid(rec):
    """
    Extract sequence ID from fasta header

    Args:
        rec (biopython record)

    Returns:
        str: sequence id that Entrez API can recognize
    """    
    
    if rec.id.startswith("UniRef100"):
        splitid = rec.id.split('_')[1]
    elif rec.id.startswith("sp|"):
        splitid = rec.id.split('|')[1]
    elif rec.id.endswith('-F1'):
        splitid = rec.id[:-3]
    elif ':(' in rec.id:
        splitid = rec.id.split(':')[0]
    else:
        splitid = rec.id
    return splitid

def store_taxids_dict(fp, fname='seqid_taxid.txt'):
    unique_taxids = {}
    with open(os.path.join(fp, fname), 'r') as f:
        content = f.readlines()
        for line in content:
            line = line.split()
            seqid = line[0]
            if not line[1] == "None":
                taxid = int(line[1])
            if not taxid in unique_taxids.keys():
                unique_taxids[taxid] = [seqid]
            else:
                unique_taxids[taxid].append(seqid)
    return unique_taxids

def records_to_fasta(records, out_fp):
    """
    Write dictionary items to fasta file

    Args:
        records (dict): dictionary of format d[sequence_id] = sequence
        out_fp (str): location where to store the fasta file
    """    
    with open(out_fp, 'w+') as out:
        for key, value in records.items():
            out.write(">" + key + "\n")
            out.write(str(value) + "\n")

def merge_dict(dict1, dict2):
    return(dict2.update(dict1))