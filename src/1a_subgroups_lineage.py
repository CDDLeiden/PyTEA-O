import os
import sys
import argparse
import pandas as pd
from collections import Counter
from obtain_lineage import filter_sequences, records_to_fasta, sequence_records

"""
Script create subgroup files, based on their lineage

    for every possible lineage split (species, family, clade, etc.) a new subfamily_{split} file is made which list the
    sequence id and the lineage
    if the sequencethreshold argument is provided, only sequenceids are kept that have a subgroup bigger than the defined threshold
    For every subfamily file a new alignment file is made, based on the sequences in the subfamily file.
"""

def lineage_fname(lin):
    return f"subfamily_{lin}.txt"

def remove_file(fp):
    if os.path.exists(fp):
        os.remove(fp)
    
def extract_from_tsf(fp, col_num):
    res = []
    with open(fp, 'r') as f:
        content = f.readlines()
        for line in content:
            line = line.split()
            res.append(line[col_num])
    return res

def length_file(fp):
    with open(fp,"r") as f:
        return len(f.readlines())

def remove_row_from_tsf(to_remove: list, fp: str):
    # Read all lines 
    with open(fp, "r") as f:
        lines = f.readlines()
    with open(fp, "w") as f:
        for line in lines:
            line = line.split()
            if not line[1] in to_remove:
                f.write(line[0] + '\t' + line[1] + "\n")

def lineage_to_file(item, df, out_fn):
    for _, row in df.iterrows():
        lineage = row[str(item)]
        if not lineage == "unknown":
            seqid = row['seqid']
            outpath = os.path.join(out_fn, lineage_fname(item))
            with open(outpath, 'a+') as file:
                file.write(str(seqid) + '\t' + str(lineage) + '\n')

def cleanup(output_location):
    if os.path.exists(os.path.join(output_location,"subfamily_Unnamed:_0.txt")):
        os.remove(os.path.join(output_location,"subfamily_Unnamed:_0.txt"))

def main(args):

    # Read lineage information
    df = pd.read_csv(args.lineage_fp)

    # Remove columns from df
    cols_to_remove = ['taxid', 'seqid']
    all_df_colums = df.columns.values.tolist()
    lineage_items = [elem for elem in all_df_colums if not elem in cols_to_remove]

    for item in lineage_items:

        # Check if output file already exists and remove before appending to old files  
        out_fp = os.path.join(args.subgroup_out_fp, lineage_fname(item))
        remove_file(out_fp)

        # Append known lineages in df to file(s). Each for every lineage split
        lineage_to_file(item, df, args.subgroup_out_fp)  

        # Count instances per subgroup and remove small subgroups below threshold
        out_fp = os.path.join(args.subgroup_out_fp, lineage_fname(item))
        subgroup_values = extract_from_tsf(out_fp, col_num=1)
        counts = Counter(subgroup_values)  # Count values and remove those below thresholds
        
        seqids_to_be_removed = []
        for value, count in counts.items():
            if count < int(args.sequencethreshold):
                seqids_to_be_removed.append(value)
        seqids_to_be_removed.append('unknown')

        remove_row_from_tsf(seqids_to_be_removed, out_fp)

        # Read sequences and extract sequence IDs to keep
        ids_to_keep = extract_from_tsf(out_fp, col_num=0)

        # Filter MSA to only contain sequences to keep
        if len(ids_to_keep) > 1 and (args.ref_id in ids_to_keep):
            all_records = sequence_records(args.msa_fp)
            filtered_sequences = filter_sequences(all_records, ids_to_keep)
            fname = os.path.basename(args.msa_fp)
            fname_new = f"{item}_{fname}"
            if not item == "Unnamed:_0":
                out_fp = os.path.join(os.path.dirname(args.msa_fp), fname_new)  # args.msa_fp
                records_to_fasta(filtered_sequences, out_fp)

        # Remove empty files
        outpath = os.path.join(args.subgroup_out_fp, lineage_fname(item))
        if length_file(outpath) == 0:
            remove_file(outpath)

    # Remove files for which reference seq does not have value
    cleanup(args.subgroup_out_fp)

def read_arguments():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("-t", "--sequencethreshold", help="minimum amount of sequences per subgroup", default=1)
    parser.add_argument("-l", "--lineage_fp", help="CSV file with lineage information")
    parser.add_argument("-s", "--subgroup_out_fp", help="folder where to store the txt files with subgroups")
    parser.add_argument("-m", "--msa_fp", help="MSA")
    parser.add_argument("-r", "--ref_id", help="reference sequence")
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    arguments = read_arguments()
    main(arguments)            