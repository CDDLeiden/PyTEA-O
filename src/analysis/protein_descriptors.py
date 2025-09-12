#!/usr/bin/env python3

import os
import sys
import numpy as np
import pandas as pd

from src.utils.general import get_file_name
from src.utils.sequence import SequenceUtilities

def load_residue_counts(summary_file:str) -> dict:

	positions_rescounts = {}

	with open(summary_file, 'r') as IN:
		
		for line in IN:
			line = line.strip()

			if line == "":
				continue

			if line[0] == "#":
				continue

			line = line.split("\t")

			msa_pos = int(line[1])
			res = line[5]

			# Count all residues for that position
			res_counts = {x.split(":")[0]:int(x.split(":")[-1]) for x in res.split(";") if x.split(":")[0]}
			# Remove gaps and 0-counts
			res_counts = {k: v for k, v in res_counts.items() if k != '-' and v > 0}
			# Store the residue counts for that position
			positions_rescounts[msa_pos] = res_counts
	
	return positions_rescounts


def calculate_descriptors(descriptor_values:dict, residue_counts:dict,outfile:str) -> None:
	"""
	Calculate the standard deviation of the descriptor values for each MSA position and scale them between 0 and 1"
	"""

	std_desc_values = {}

	for msa_pos, res_counts in residue_counts.items():

		# Initialize dictionary for every MSA position
		tmp = {key: np.array([]) for key in descriptor_values['A'].keys()}

		for res, count in res_counts.items():
			if res not in descriptor_values:
				print(f"Residue {res} not found in descriptor values. Skipping...")
				continue

			# Get the descriptor values for that residue
			descriptor_values_res = descriptor_values[res]
			for desc_key, desc_value in descriptor_values_res.items():
				tmp[desc_key] = np.append(tmp[desc_key], np.full(count, desc_value))

		# Calculate the std for each descriptor
		std_desc_values[msa_pos] = {key: round(np.std(value),3) for key, value in tmp.items()}

	# Convert to DataFrame
	std_desc_values = pd.DataFrame.from_dict(std_desc_values, orient='index')
	# Scale descriptor values between 0 and 1
	std_desc_values_scaled = ((std_desc_values - std_desc_values.min()) / (std_desc_values.max() - std_desc_values.min())).round(3)
	# Save to TSV file
	desc_names = std_desc_values.keys().to_list()
	with open(f"{outfile}",'w') as OUT:
		OUT.write(f"## MSA_Pos\t" + "\t".join(desc_names) + "\n")
		for _, row in std_desc_values_scaled.iterrows():
			OUT.write(f"{int(row.name)}\t" + "\t".join([str(row[desc]) for desc in desc_names]) + "\n")


def run(descriptor_file:str,summary_file:str,outdir:str) -> None:
	
	if not os.path.exists(outdir):
		os.makedirs(outdir)

	descriptors = SequenceUtilities.Sandberg_Zscales
	res_counts = load_residue_counts(summary_file=summary_file)
	exit()

	desc_name = get_file_name(descriptor_file)
	outfile = f"{outdir}/prot_descriptors_{desc_name}.tsv"

	calculate_descriptors(descriptor_values=descriptors, residue_counts=res_counts, outfile=outfile)

if __name__ == "__main__":

	from argparse import ArgumentParser
	from os.path import dirname,abspath

	parser = ArgumentParser()

	parser.add_argument("-o","--outdir",required=False,type=str,default="zscales")
	parser.add_argument("-d","--descriptor_file",required=False,default=f"{dirname(abspath(__file__))}/../DB_files/ZscaleSandberg.txt")
	parser.add_argument("-s","--summary_file",required=True,default=None,type=str,help="shannon_entropy.summary file")

	args = parser.parse_known_args()[0]
	
	run(descriptor_file=args.descriptor_file, 
		summary_file=args.summary_file, 
		outdir=args.outdir)