#!/usr/bin/env python3

import os
import sys
import numpy as np
import pandas as pd
from prodec import ProteinDescriptors
from convert_alignment import read_MSA


def calculate_protein_descriptor(aligned_sequences:list, descriptor:str) -> pd.DataFrame:
	"""
	Calculate protein descriptors based on a list of aligned sequences
	"""
	desc_factory = ProteinDescriptors()
	print(f"Calculating {descriptor} descriptor")
	prodec_descriptor = desc_factory.get_descriptor(descriptor)
	# result = prodec_descriptor.pandas_get(aligned_sequences, gaps='omit')
	result = prodec_descriptor.pandas_get(aligned_sequences, gaps=0.0)
	return result

def scale_df(arr):
	df = pd.DataFrame(arr)
	df_scaled = (df - df.min()) / (df.max() - df.min())
	return df_scaled

def std(df):
	return df.std(axis=0)

def reshape_df(df):
	arr = np.array(df)
	matrix = np.zeros((5, int(len(df)/5)))  # 5 different Z-scales (Sandberg)
	count = 0
	for y in range(int(len(df)/5)):
		for x in range(5):
			matrix[x][y] = arr[count]
			count += 1
	return matrix

def run(msa_file:str,outdir:str="zscales",descriptors:str='Zscale Sandberg'):


	AAs = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']

	# Create output directory
	if not os.path.exists(outdir):
		os.makedirs(outdir)
	
	prefix,msa = read_MSA(msa_file)

	sequences = ["".join([res if res.upper() in AAs else '-' for res in sequence]) for sequence in msa.values()]

	# Calculate Z-scales
	z_scales = calculate_protein_descriptor(sequences, 'Zscale Sandberg')

	# Replace the 0s (gaps) by NaN
	z_scales.replace(0, np.nan, inplace=True)
	# Calculate std for each z-scale
	z_scales_std = std(z_scales)
	# Reshape
	z_scales_reshaped = reshape_df(z_scales_std)
	# Scale z-scales
	z_scales_scaled = scale_df(z_scales_reshaped.T)
	# Round to 2 decimals
	z_scales_scaled = z_scales_scaled.round(2)
	z_scales_scaled.rename(columns={0: 'Z1', 
									1: 'Z2',
									2: 'Z3',
									3: 'Z4',
									4: 'Z5'}, inplace=True)
	
	# Add msa positions to dataframe
	msa_positions = [i for i,_ in enumerate(sequences[0])]
	z_scales_scaled = z_scales_scaled.apply(pd.to_numeric, errors='coerce')
	z_scales_scaled['MSA_position'] = msa_positions

	# Save Zscales to file
	with open(f"{outdir}/zscales.tsv",'w') as OUT:
		OUT.write(f"## MSA_position\tZ1\tZ2\tZ3\tZ4\tZ5\n")
		for _, row in z_scales_scaled.iterrows():
			OUT.write(f"{int(row['MSA_position'])}\t{row['Z1']}\t{row['Z2']}\t{row['Z3']}\t{row['Z4']}\t{row['Z5']}\n")


if __name__ == "__main__":

	from argparse import ArgumentParser

	GetOptions = ArgumentParser()

	GetOptions.add_argument("-m","--msa_file",required=True,help="MSA in .fasta format")
	GetOptions.add_argument("-o","--outdir",required=False,type=str,default="zscales")
	GetOptions.add_argument("-d","--descriptors",required=False,default='Zscale Sandberg')

	args = GetOptions.parse_known_args()[0]

	run(msa_file=args.msa_file,outdir=args.outdir,descriptors=args.descriptors)