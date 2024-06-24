import os
import sys
import time
import argparse
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio import Entrez

"""
Script to obtain lineage of sequences in MSA using the Entrez API.

A new MSA is created as output that only stores the sequences for which the lineage could be found.
The lineage information is stored in a CSV file with _lineage.csv extension
"""

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
        handle.close()
        tax = int(taxonList[0]["TaxId"])
        if tax == 1:
            print(code, tax)
    except RuntimeError:
        tax = None
    return tax

def Entrez_lineage_from_taxid(tax_num):
    """
    Obtain lineage from Entrez and return dictionary with format e.g.
    
        {500: [('no rank', 'cellular organisms'),
               ('superkingdom', 'Bacteria'),
               ('clade', 'Terrabacteria group'),
               ('phylum', 'Chloroflexi'),
               ('class', 'Thermomicrobia'),
               ('order', 'Thermomicrobiales'),
               ('family', 'Thermomicrobiaceae'),
               ('genus', 'Thermomicrobium')],
    """
    tax_dict = {}
    tax_dict[tax_num] = []
    handle = Entrez.efetch(db="taxonomy", id=tax_num, retmode="xml")
    time.sleep(0.33)
    recordList = Entrez.read(handle)
    handle.close()
    for record in recordList:
        if record["Lineage"]:
            res = record["LineageEx"]
            for i in res:
                for key, value in i.items():
                    if key == "ScientificName":
                        name = value
                    elif key == "Rank":
                        rank = value
                tax_dict[tax_num].append((rank, name))
    return tax_dict
    
def format_lineage(entrez_lineage_info):
    """
    Format Entrez lineage information and store in pd.DataFrame
    """    
    # Order of lineage:
    # Superkingdom, Clade, Phylum, Class, Order, Family, Genus, Species
    
    superkingdom = list()
    clade = list()
    phylum = list()
    Class = list()
    order = list()
    family = list()
    genus = list()
    species = list()
    seqid = list()
    taxid = list()

    for key, entrez_lineage in entrez_lineage_info.items():

        t_superkingdom = []
        t_clade = []
        t_clade = []
        t_phylum = []
        t_class = []
        t_order = []
        t_family = []
        t_genus = []
        t_species = []

        try:
            seqid.append(key[0])
            taxid.append(key[1])
        except TypeError:  # Just a single key
            taxid.append(key)

        for val in entrez_lineage:
            if val[0] == 'superkingdom':
                t_superkingdom.append(val[1])
            elif val[0] == 'phylum':
                t_phylum.append(val[1])
            elif val[0] == 'class':
                t_class.append(val[1])
            elif val[0] == 'order':
                t_order.append(val[1])
            elif val[0] == 'family':
                t_family.append(val[1])
            elif val[0] == 'genus':
                t_genus.append(val[1])
            elif val[0] == 'species':
                t_species.append(val[1])
            elif val[0] == 'clade':
                t_clade.append(val[1])

        # If taxonomy is not obtained add 'unknown'
        if len(t_superkingdom) == 0:
            superkingdom.append('unknown')
        if len(t_phylum) == 0:
            phylum.append('unknown')
        if len(t_clade) == 0:
            clade.append('unknown')
        if len(t_class) == 0:
            Class.append('unknown')
        if len(t_order) == 0:
            order.append('unknown')
        if len(t_family) == 0:
            family.append('unknown')
        if len(t_genus) == 0:
            genus.append('unknown')
        if len(t_species) == 0:            
            species.append('unknown')

        if len(t_superkingdom) == 1:
            superkingdom.append(t_superkingdom[0])
        if len(t_phylum) == 1:
            phylum.append(t_phylum[0])
        if len(t_clade) == 1:
            clade.append(t_clade[0])
        if len(t_class) == 1:
            Class.append(t_class[0])
        if len(t_order) == 1:
            order.append(t_order[0])
        if len(t_family) == 1:
            family.append(t_family[0])
        if len(t_genus) == 1:
            genus.append(t_genus[0])
        if len(t_species) == 1:            
            species.append(t_species[0])

        # If there are multiple options per lineage_part > merge
        if len(t_superkingdom) > 1:
            merged = '_'.join(t_superkingdom)
            superkingdom.append(merged)
        if len(t_phylum) > 1:
            merged = '_'.join(t_phylum)
            phylum.append(merged)
        if len(t_clade) > 1:
            merged = '_'.join(t_clade)
            clade.append(merged)
        if len(t_class) > 1:
            merged = '_'.join(t_class)
            Class.append(merged)
        if len(t_order) > 1:
            merged = '_'.join(t_order)
            order.append(merged)
        if len(t_family) > 1:
            merged = '_'.join(t_family)
            family.append(merged)
        if len(t_genus) > 1:
            merged = '_'.join(t_genus)
            genus.append(merged)
        if len(t_species) > 1:
            merged = '_'.join(t_species)
            species.append(merged)  
            
    df = pd.DataFrame()
    try:
        df['seqid'] = seqid
    except:
        df['seqid'] = None
    df['taxid'] = taxid
    df['superkingdom'] = superkingdom
    df['phylum'] = phylum
    df['clade'] = clade
    df['class'] = Class
    df['order'] = order
    df['family'] = family
    df['genus'] = genus
    df['species'] = species
    return df

def get_output_filepath(fp, extension):
    fname = os.path.basename(fp)
    fname_no_ext = os.path.splitext(fname)[0]
    out_fp = fname_no_ext + extension
    return out_fp

def filter_sequences(all_records, records_to_keep):
    """
    Keep subset of sequences based on list with ids to keep

    Args:
        all_records (biopython records): _description_
        records_to_keep (list): list of IDs to keep

    Returns:
        dict: dictionary of format d[sequence_id] = sequence
    """    
    filtered_sequences = {}
    for record in all_records:
        seqid = seqid_from_fullid(record)
        if seqid in records_to_keep:
            filtered_sequences[seqid] = record.seq
    return filtered_sequences

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

def sequence_records(msafile, extension="fasta"):
    records = list(SeqIO.parse(msafile, extension))
    return records

def main(args):

    records = sequence_records(args.msafile)

    taxonomy_combined = {}
    for record in records:

        seqid = seqid_from_fullid(record)  # Extract seqid that Entrez can recognize
        tax_num = run_entrez_api(seqid)  # Get taxonomy number
        print(seqid, tax_num)
        if tax_num:
            try:
                result = Entrez_lineage_from_taxid(tax_num)
                taxonomy_combined[seqid, tax_num] = result[tax_num]
            except:
                print(f"Could not obtain lineage for {seqid}")
    
    # Store lineage information in pd.DataFrame and save
    df = format_lineage(taxonomy_combined)
    out_fname = get_output_filepath(args.msafile, "_lineage.csv")
    df.to_csv(os.path.join(args.output_location, out_fname), index=False)

    # Make new MSA which only contains sequences for which data could be obtained from entrez
    seqids_with_taxonomy = [key[0] for key, _ in taxonomy_combined.items()]
    filtered_sequences = filter_sequences(records, seqids_with_taxonomy)
    out_fname = get_output_filepath(args.msafile, "_filtered.fasta")
    out_fp = os.path.join(args.output_location, out_fname)
    records_to_fasta(filtered_sequences, out_fp)
    
def read_arguments():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("-i", "--msafile", help="Multiple Sequence Alignment in FASTA format")
    parser.add_argument("-o", "--output_location", help="file location where output should be stored")
    parser.add_argument("-e", "--email", help="Email adress that is used by Entrez API")
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    arguments = read_arguments()
    main(arguments)    