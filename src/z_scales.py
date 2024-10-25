#!/usr/bin/env python3

import os
import sys
import numpy as np
import pandas as pd
import seaborn as sns
from Bio import AlignIO
import matplotlib.pyplot as plt
from prodec import ProteinDescriptors
from sklearn.preprocessing import MinMaxScaler


def calculate_protein_descriptor(aligned_sequences:list, descriptor:str) -> pd.DataFrame:
	"""
	Calculate protein descriptors based on a list of aligned sequences
	"""
	desc_factory = ProteinDescriptors()
	prodec_descriptor = desc_factory.get_descriptor(descriptor)
	result = prodec_descriptor.pandas_get(aligned_sequences)
	return result

def read_alignment(fp, extension="fasta"):
    return AlignIO.read(fp, extension)

def scale_df(arr):
	df = pd.DataFrame(arr)
	df.head()
	scaler = MinMaxScaler()
	df_scaled = pd.DataFrame(scaler.fit_transform(df), columns=df.columns)
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

# def heatmap(df_zscales, xlabels):
#     fig, ax = plt.subplots(figsize=(5, 5))
#     sns.heatmap(df_zscales, ax=ax, cmap='Reds', annot=False, cbar=False)
#     ax.xaxis.tick_bottom()
#     ax.patch.set_alpha(0)
#     ax.set_xlabel('Residue')
#     ax.set_ylabel('Z-scale')
#     # put the major ticks at the middle of each cell
#     ax.set_xticks(np.arange(df_zscales.shape[1]) + 0.5, minor=False)
#     ax.set_yticks(np.arange(df_zscales.shape[0]) + 0.5, minor=False)
#     return fig, ax

def get_residues_per_position(msa):
	positions = {}
	for num, record in enumerate(msa):
		positions[num] = [aa for aa in record.seq]
	return positions

def run(args=None):

	# Create output directory
	if not os.path.exists(args.outdir):
		os.makedirs(args.outdir)
	
	# desc_factory = ProteinDescriptors()
	msa = read_alignment(args.msa_file)
	# Get positions from the msa
	positions = get_residues_per_position(msa)
	# Calculate Z-scales
	truncated_sequences = []
	for _, v in positions.items():
		truncated_sequences.append(''.join(v))
	z_scales = calculate_protein_descriptor(truncated_sequences, 'Zscale Sandberg')
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
	msa_positions = [i for i in range(len(msa[0].seq))]
	z_scales_scaled = z_scales_scaled.apply(pd.to_numeric, errors='coerce')
	z_scales_scaled['MSA_position'] = msa_positions

	# Save Zscales to file
	z_scales_scaled.to_csv(f'{args.outdir}/zscales.txt', index=False, sep='\t')

	if args.shannon_entropy_file:
		# Read shannon entropy file
		shannon_entropy = pd.read_csv(args.shannon_entropy_file, sep='\t')
		se_msa_positions = shannon_entropy['MSA_Pos'].tolist()
		se_res = shannon_entropy['Residue'].tolist()
		se_seq_pos = shannon_entropy['## Sequence_Pos'].tolist()
		# Only keep the positions that are in the shannon entropy file
		z_scales_scaled = z_scales_scaled[z_scales_scaled['MSA_position'].isin(se_msa_positions)]
		z_scales_scaled['Residue'] = se_res
		z_scales_scaled['Sequence_Pos'] = se_seq_pos
		seqid = args.shannon_entropy_file.split('/')[-1].split('.')[0]
		z_scales_scaled.to_csv(f'{args.outdir}/{seqid}.zscales.txt', index=False, sep='\t')

if __name__ == "__main__":

	from argparse import ArgumentParser

	GetOptions = ArgumentParser()

	GetOptions.add_argument("-m","--msa_file",required=True)
	GetOptions.add_argument("-o","--outdir",required=False,type=str,default="zscales")
	GetOptions.add_argument("-d","--descriptors",required=False,default='Zscale Sandberg')
	GetOptions.add_argument("-s","--shannon_entropy_file",required=False,default=None)
	# GetOptions.add_argument("-r","--reference_id",required=False,default=None)

	run(GetOptions.parse_known_args()[0])