#!/usr/bin/env python3

import os
import sys
from math import log
from statistics import median
import argparse
import numpy as np
import pandas as pd
from time import time
from typing import Tuple
from src.grouping.upgma_subfam_grouping import run as run_upgma
from src.utils.msa import get_reference_indices, read_MSA, get_residue_counts, load_msa
from src.utils.general import read_subfamilies,create_directory
from src.utils.multiprocess import get_multiprocessing_module
from src.io.output import write_references,write_shannon_entropy_temp_results

def calculate_shannon_entropy(msa:pd.DataFrame,weigh:bool=True) -> Tuple[np.ndarray,np.ndarray,np.ndarray]:

	"""
	Calculates the weighted Shannon Entropy for the provided MSA utilizing the proposed equation by Ye, 2008

	Parameters
	msa (pd.DataFrame): DataaFrame containing all residue counts for all sequence indexes for all sequences in the MSA

	Output
	E_i (np.ndarray): Weighted Shannon Entropy for all sequence indexes in the MSA
	a_i (np.ndarray): Residue counts for all sequence indexes in the MSA
	G_i (np.ndarray): Gap counts for all sequence indexes in the MSA
	"""

	n_j,_ = msa.shape

	_,num_of_MSA_res = msa.loc[msa.index[0],'array'].shape

	## Occurrences of residues at position of i
	residue_counts:dict = get_residue_counts(msa=msa)

	## Occurrences of AA at position i in the subfamily, minus gaps, normalized by number of proteins in subfamily
	## Array is of size number_of_residues x 20 (number of amino acids)
	a_i:np.array = residue_counts[0:-1]/n_j

	## ln(0) ==> bad things... instead make ln(1) ==> 0
	a_i[a_i==0] = 1

	## Number of gaps at position in the subfamily
	G_i:np.array = residue_counts[-1]

	G_i = G_i.reshape(1,num_of_MSA_res)

	## Normalization factor
	W_j:float = 1/np.log(n_j) if n_j < 20 else 1/np.log(20)

	W_j = W_j if weigh else 1/np.log(n_j)

	## Add weighted subfamily entropy to total entropy
	E_i = W_j*((-np.sum(a_i*np.log(a_i),axis=0)).reshape(1,num_of_MSA_res)-G_i*(1/n_j)*np.log(1/n_j))

	return E_i,residue_counts[0:-1],G_i

def superfamily_shannon_entropy(msa:pd.DataFrame,outdir:str,file_prefix:str="TEA") -> Tuple[np.ndarray,np.ndarray,np.ndarray]:

	outfile = f"{outdir}/{file_prefix}.superfamily_entropy"

	_,num_of_MSA_res = msa.loc[msa.index[0],'array'].shape

	E_i,a_i,G_i = calculate_shannon_entropy(msa=msa,weigh=False)

	with open(outfile,'w') as OUT:
		OUT.write(f"## MSA_position\tShannon_Entropy\tSequences_at_pos\tNum_of_gaps\n")
		for index in range(num_of_MSA_res):

			sh_entropy = E_i[0][index]
			non_gapped_res = int(np.sum(a_i[:,index]))
			num_of_gaps = int(G_i[0][index])
			
			# MSA Position
			OUT.write(f"{index}")

			## Shannon Entropy
			OUT.write(f"\t{sh_entropy: >.2f}")

			## Sequences at position
			OUT.write(f"\t{non_gapped_res}")

			## Number of gaps
			OUT.write(f"\t{num_of_gaps}\n")

	return E_i,a_i,G_i

def shannon_entropy_multiprocess(msa:pd.DataFrame,temp_entropy_file:str,results:dict,subfamilies:dict,keys:list,threads:int=1):
	
	def worker(key,subfamily_groupings,queue,msa):
		entropy = calculate_subfamily_shannon_entropy(msa,subfamily_groupings)
		queue.put([key,entropy])
		return entropy

	mp = get_multiprocessing_module()

	queue = mp.Manager().Queue()
	write_temp_process = mp.Process(target=write_shannon_entropy_temp_results,args=(queue,temp_entropy_file))
	write_temp_process.start()

	with mp.Pool(processes=threads) as executor:
		for subfamily in keys:
			if subfamily in results.keys():
				continue
			executor.apply_async(worker,args=(subfamily,subfamilies[subfamily],queue,msa))
		executor.close()
		executor.join()

	queue.put(False)
	write_temp_process.join()

def calculate_subfamily_shannon_entropy(msa:pd.DataFrame,groupings:dict) -> Tuple[int,np.ndarray]:

	"""Returns the average entropy of all i-th positions within an MSA according to Ye's subfamily weighted averaging.

	### Input:
		**groupings_data**: _tuple_ (branch_index,subfamily_groupings)
		Subfamily groupings, where a key represents the family ID and the value is an array of sequence IDs from the MSA.
	"""

	_,num_of_MSA_res = msa.loc[msa.index[0],'array'].shape

	se = np.zeros((1,num_of_MSA_res))

	start = time()

	for subfamily_group_number in groupings.keys():

		## Number of proteins in the subfamily
		n_j:int = len(groupings[subfamily_group_number])

		## A subfamily with a single sequence does not provide information on the Two-Entropy approach
		if n_j == 1:
			continue

		E_i,_,_ = calculate_shannon_entropy(msa=msa.loc[groupings[subfamily_group_number]])

		se += E_i

	## Return residue-wise entropy normalized by number of subgroups
	return se/len(groupings.keys())