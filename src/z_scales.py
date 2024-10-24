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


def calculate_protein_descriptor(aligned_sequences, descriptor):
	"""
	Calculate protein descriptors of choice for unique sequences in the bioactivity dataset

	Parameters
	----------
	aligned_sequences : list
		List of aligned sequences
	descriptor : str
		Descriptor to calculate

	Returns
	-------
	pandas.DataFrame
		Dataset with sequence and features for the protein descriptors of interest for sequences in the bioactivity dataset
	"""
	desc_factory = ProteinDescriptors()
	prodec_descriptor = desc_factory.get_descriptor(descriptor)
	result = prodec_descriptor.pandas_get(aligned_sequences)
	return result


def truncate_seq(msa, pos):
	"""
	Truncate sequences in MSA to only include residues at positions of interest

	Parameters
	----------
	msa : Bio.Align.MultipleSeqAlignment
		Multiple sequence alignment
	pos : list
		List of positions of interest

	Returns
	-------
	list
		List of truncated sequences
	"""
	sequences = {}
	for record in msa:
		sequences[record.id] = []

	for p in pos:
		for record in msa:
			sequences[record.id].append(record.seq[p])

	truncated_seqs = []
	recids = []
	for k, v in sequences.items():
		sequences[k] = ''.join(v)
		truncated_seqs.append(sequences[k])
		recids.append(k)

	check_seq(truncated_seqs)  # Check if sequences contain only valid amino acids
	
	return truncated_seqs, recids


def check_seq(seq, print_invalid=False):
	"""
	Check if sequence contains only valid amino acids

	Parameters
	----------
	seq : list of sequences

	Returns
	-------
	list
	"""
	aas = validaas()
	incorrect_positions = []
	for num, i in enumerate(seq):
		if not all(aa in aas for aa in i):
			if print_invalid:
				print('Invalid amino acid sequence found. Please check your input file.')
			incorrect_positions.append(num)
			# Find out which sequence is invalid
			for j in i:
				if j not in aas:
					if print_invalid:
						print('Invalid amino acid:', j)
						print('Invalid sequence:', i)
						print('Invalid sequence number:', num)
						print('Please check your input file.')

	return incorrect_positions

def impute_seq(seqlist, type='remove'):
	"""
	Impute missing residues in MSA

	Parameters
	----------
	sequences : list of sequences
	pos : int
	type : str (default: 'remove') Remove gaps and replace with NaN value

	Returns
	-------
	list
		List of imputed sequences
	"""
	# TODO FIX THIS
	columns = len(seqlist[0])
	colnames = [f'pos_{i}' for i in range(columns)]
	# TODO Add some statistics about what residues are imputed and what is the mode
	df_seqlist = pd.DataFrame(seqlist, columns=['seq'])
	df = df_seqlist['seq'].apply(lambda x: pd.Series(list(x)))
	df.replace('-', np.nan, inplace=True)
	return df
	# df['seq'] = df[:].agg(''.join, axis=1)
	# imputed_seqlist = df['seq'].tolist()
	# return imputed_seqlist


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
	matrix = np.zeros((5, int(len(df)/5)))
	count = 0
	for y in range(int(len(df)/5)):
		for x in range(5):
			matrix[x][y] = arr[count]
			count += 1
	return matrix  # 5 different Z-scales (Sandberg)


def validaas():
	return ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'W', 'Y', 'V', 'M', 'N', 'P', 'Q', 'R', 'S', 'T']


def heatmap(df_zscales, xlabels):

    fig, ax = plt.subplots(figsize=(5, 5))
    sns.heatmap(df_zscales, ax=ax, cmap='Reds', annot=False, cbar=False)
    ax.xaxis.tick_bottom()
    ax.patch.set_alpha(0)
    ax.set_xlabel('Residue')
    ax.set_ylabel('Z-scale')

    # put the major ticks at the middle of each cell
    ax.set_xticks(np.arange(df_zscales.shape[1]) + 0.5, minor=False)
    ax.set_yticks(np.arange(df_zscales.shape[0]) + 0.5, minor=False)

    return fig, ax

def read_msa(msa_file):
	msa = AlignIO.read(msa_file, 'fasta')
	return msa

def get_residues_per_position(msa):
	positions = {}
	for num, record in enumerate(msa):
		positions[num] = [aa for aa in record.seq]
	# for i in range(len(msa[0].seq)):
	# 	positions[i] = [record.seq[i] for record in msa]
	return positions

def get_mapped_positions(msa):
	positions = {}
	for num in range(len(msa[0].seq)):
		positions[num] = [record.seq[num] for record in msa]
	return positions

def remove_gaps(positions):
	# Remove gaps from positions
	positions = {k: [aa.replace('-', '') for aa in v] for k, v in positions.items()}
	# Remove positions with invalid amino acids
	# TODO Check this
	# positions = {k: ''.join(aa for aa in v if aa in validaas()) for k, v in positions.items()}
	return positions


def run(args=None):

	#TODO Make output dir if it does not exist
	if not os.path.exists(args.outdir):
		os.makedirs(args.outdir)
	
	desc_factory = ProteinDescriptors()
	# print('Available ProDEC descriptors: ', desc_factory.available_descriptors)  # Print available descriptors
	msa = read_msa(args.msa_file)
	
	# Get positions from the msa
	positions = get_residues_per_position(msa)

	mapped_positions = get_mapped_positions(msa)
	mapped_cleaned_positions = remove_gaps(mapped_positions)

	for key, value in mapped_positions.items():
		print(key, len(value), value)

	print('######################################################')

	for key, value in mapped_cleaned_positions.items():
		result = value.count('')
		print(key, len(value) - result)

	# Remove gaps and unvalid amino acids
	# positions = remove_gaps(positions)

	# for k, v in positions.items():
	# 	print(k, v)

	# Calculate Z-scales
	truncated_sequences = []
	positions2 = {}
	for k, v in positions.items():
		truncated_sequences.append(''.join(v))
		positions2[k] = ''.join(v)

	# for num, i in enumerate(truncated_sequences):
	# 	print(num, len(i), i)

	print("truncated_sequences", len(truncated_sequences), truncated_sequences[0])

	z_scales = calculate_protein_descriptor(truncated_sequences, 'Zscale Sandberg')
	print("z_scales_sandberg", z_scales.shape)

	# Replace the 0s (gaps) by NaN
	z_scales.replace(0, np.nan, inplace=True)

	# Calculate std for each z-scale
	z_scales_std = std(z_scales)
	print("z_scales_std", z_scales_std.head())

	# Reshape and save
	z_scales_reshaped = reshape_df(z_scales_std)
	print("z_scales_reshaped", z_scales_reshaped.shape)

	# Scale z-scales
	z_scales_scaled = scale_df(z_scales_reshaped.T)
	print("z_scales_scaled", z_scales_scaled.shape)

	# Round to 2 decimals
	z_scales_scaled = z_scales_scaled.round(2)

	z_scales_scaled.rename(columns={0: 'Z1', 
									1: 'Z2',
									2: 'Z3',
									3: 'Z4',
									4: 'Z5'}, inplace=True)
	
	# TODO Add residues to dataframe (if highlight_residues is provided) otherwise add position numbers
	positions = [i for i in range(len(msa[0].seq))]
	# print("positions", positions)

	# Transpose
	z_scales_scaled = z_scales_scaled.apply(pd.to_numeric, errors='coerce')
	print("z_scales_scaled", z_scales_scaled.shape)
	z_scales_scaled['MSA_position'] = positions
	# z_scales_scaled = z_scales_scaled.iloc[::-1] # Reverse order so Z1 is at the bottom

	# TODO Add residues
	if args.highlight_residues:
		for rec in msa:
			if rec.id == args.highlight_residues:
				highlighted_residues = rec.seq
				print("highlighted_residues", highlighted_residues)
				z_scales_scaled['Res'] = [aa for aa in highlighted_residues]

	# Save Zscales to CSV
	z_scales_scaled.to_csv(f'{args.outdir}/Zscales.csv', index=False, sep=';')

if __name__ == "__main__":

	from argparse import ArgumentParser

	GetOptions = ArgumentParser()

	GetOptions.add_argument("-s","--shannon_entropy_file",required=True,type=str)
	GetOptions.add_argument("-o","--outdir",required=False,type=str,default="zscales")
	GetOptions.add_argument("-y","--highlight_residues",required=False,type=str)
	GetOptions.add_argument("-c","--consensus_sequence_file",required=False,type=str)
	GetOptions.add_argument("-a","--average_entropy_file",required=False)
	GetOptions.add_argument("-d","--descriptors",required=False,default='Zscale Sandberg')
	GetOptions.add_argument("-m","--msa_file",required=True)

	run(GetOptions.parse_known_args()[0])