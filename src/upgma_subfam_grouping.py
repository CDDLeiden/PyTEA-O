#!/usr/bin/env python3

from os.path import isfile
from os import makedirs
import sys
import argparse
import numpy as np
import pandas as pd
from time import time
from convert_alignment import read_MSA
from warnings import simplefilter
import gc

simplefilter(action='ignore',category=pd.errors.PerformanceWarning)

## MacOS (darwin) utilizes a different multiprocessing library than Linux and Windows. Need to check which OS
## OS is running to prevent SNAFUs
multiprocessing = None
if 'linux' in sys.platform:
	multiprocessing = 'multiprocessing'
elif 'darwin' in sys.platform:
	multiprocessing = 'multiprocess'
elif 'win' in sys.platform:
	multiprocessing = 'multiprocessing'
mp = __import__(multiprocessing)

def generate_distance_matrix(msa:dict) -> pd.DataFrame:

	"""Returns an NxN numpy matrix with values set to -1.

	Inputs:
		msa [dict]: Dictionary of MSAs, where keys are the loci and values are the sequences

	Returns:
		np.array: NxN array with each element being set to -1
	"""

	keys:dict.keys = sorted(msa.keys())

	return pd.DataFrame(np.zeros([len(msa),len(msa)]),columns=keys,index=keys)

def distance(sequence_a_loc):

	sequence_a = sequences[sequence_a_loc]

	distances = {x:0 for x in sequences.keys()}

	for sequence_b_loc in sorted(sequences.keys()):
		if sequence_b_loc == sequence_a_loc:
			distances[sequence_b_loc] = np.inf
			break
		distances[sequence_b_loc] = sum(1 for a,b in zip(sequence_a,sequences[sequence_b_loc]) if a!=b)

	return [sequence_a_loc,distances]

def calculate_sequence_difference(msa:dict,outfile:str) -> pd.DataFrame:

	"""Calculates the distances between sequences in an MSA

	Returns:
		pandas.Dataframe: NxN lower triangular matrix
	"""

	if isfile(outfile):

		print(f"\n\t\t[N] Loading previously calculated pairwise distances...\n")
		return pd.read_csv(outfile,header=0,index_col=0)

	distance_matrix = generate_distance_matrix(msa)

	keys:dict.keys = sorted(msa.keys())
	global sequences
	sequences = msa
	
	start = time()

	results:list
	with mp.Pool(mp.cpu_count()) as executor:
		results = executor.map(distance,keys)

	print(f"\n\t\tDistance calculations took {int(time()-start)}s\n")

	pause = time()
	item:int = 0
	for sequence_loc_a,result in results:
		item += 1
		distance_matrix.loc[sequence_loc_a] = result
	print(f"\n\t\tDistance assignment took {int(time()-pause)}s\n")

	print(f"\n\t\tGenerated {len(keys)}x{len(keys)} distance matrix in {int(time()-start)}s\n")

	distance_matrix.to_csv(outfile)

	return distance_matrix

def UPGMA(distance_matrix:pd.DataFrame) -> dict:

	print(f"\n\t\tStarting UPGMA tree building")

	start = time()

	distance_matrix = distance_matrix + distance_matrix.T

	distances:dict = {x:0 for x in distance_matrix.index}

	tree:dict = {0:distance_matrix.columns.to_list()}

	while(len(distance_matrix)>2):

		flattend_index =  np.argmin(distance_matrix)

		row_i = flattend_index//distance_matrix.shape[1]
		col_i = flattend_index%distance_matrix.shape[1]

		row = distance_matrix.columns[row_i]
		col = distance_matrix.columns[col_i]

		pairwise_dist = distance_matrix.loc[row,col]

		data_A = np.delete(distance_matrix.loc[row].to_numpy(),[row_i,col_i])
		weight_A = len(row.split(","))
		data_B = np.delete(distance_matrix.loc[col].to_numpy(),[row_i,col_i])
		weight_B = len(col.split(","))

		new_data = np.append((((data_A*weight_A) + (data_B*weight_B))/(weight_A+weight_B)),np.inf)

		distance_matrix = distance_matrix.drop(index=[row,col],columns=[row,col])

		joined_name = f"{row},{col}"

		tree[pairwise_dist/2] = distance_matrix.columns.to_list() + [joined_name]

		distance_matrix[joined_name] = new_data[:-1]
		distance_matrix.loc[joined_name] = new_data

		gc.collect()

	row,col = distance_matrix.columns.to_numpy()

	tree[distance_matrix.loc[row,col]] = [",".join(tree[0])]
	
	print(f"\n\t\t\tCondensed {len(tree[0])} sequences in {int(time()-start)}s\n")

	return tree

def generate_subgroups(tree:dict,outfile:str)-> None:
	
	subfamilies = {}
	for branchpoint, distance in enumerate(sorted(tree.keys())):
		subfamilies[branchpoint] = {}
		for group_num,group in enumerate(sorted(tree[distance])):
			for member in group.split(","):
				subfamilies[branchpoint][member] = hex(group_num)

	num_of_groups = np.inf
	with open(outfile, "w") as file:
		headers = ";".join(tree[0])
		file.write(f"{headers}\n")
		for branch in subfamilies:
			grouping_line = ";".join([subfamilies[branch][member] for member in sorted(subfamilies[branch].keys())])
			file.write(f"{grouping_line}\n")
			num = len(np.unique([subfamilies[branch][member] for member in sorted(subfamilies[branch].keys())]))
			if num > num_of_groups: print("OH NO!!!")
			num_of_groups = num

	return None

def pass_subfamilies_file_check(outfile:str,labels:list) -> bool:

	if not isfile(outfile):

		return False

	subfamilies_data:pd.DataFrame

	try:

		subfamilies_data = pd.read_csv(outfile,header=0,index_col=0,delimiter=";")

	except pd.errors.EmptyDataError as error:

		print(f"\n\t\t[N] {outfile} exists, but contains no subfamily information.\n")

		return False

	if (labels != sorted(subfamilies_data.columns.to_list())) and (np.unique(subfamilies_data.iloc[-1].to_numpy()) != ['0x0']):

		print(f"\n\t\t[N] {outfile} exists, but subfamily information is incomplete.\n")

		return False

	print(f"\n\t\t[N] {outfile} exists, and subfamily information is complete.\n")
	print(f"\n\t\t[N] Terminating UPGMA grouping. Use -f|--force to force algorithm run.\n")

	return True

def run(msa_file:str,outdir:str) -> None:

	start = time()

	makedirs(outdir,0o700,exist_ok=True)

	## Read in MSA
	file_prefix,msa = read_MSA(msa_file=msa_file)

	## Create and fill initial distance matrix
	distance_matrix:pd.DataFrame = calculate_sequence_difference(msa=msa,outfile=f"{outdir}/{file_prefix}.dist")

	subfamilies_file = f"{outdir}/{file_prefix}.subfamilies"

	## No re-runs
	if pass_subfamilies_file_check(outfile=subfamilies_file,labels=sorted(msa.keys())): return None

	## Perform UPGMA tree analysis
	tree = UPGMA(distance_matrix=distance_matrix)

	## Identify subgroups from UPGMA collapse pattern
	generate_subgroups(tree=tree,outfile=subfamilies_file)

	return None

if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	parser.add_argument("-m", "--msa_file", required=True)
	parser.add_argument("-o","--outdir",default="UPGMA_SUBFAMILY_GROUPINGS")
	args = parser.parse_args()

	run(msa_file=args.msa_file,outdir=args.outdir)