#!/usr/bin/env python3

import os
import sys
from math import log
from statistics import median
import argparse
import numpy as np
import pandas as pd
from os import makedirs
from os.path import exists,isfile
from time import time
from convert_alignment import read_MSA
from typing import Tuple
from upgma_subfam_grouping import run as run_upgma

multiprocessing = None
if 'linux' or 'win' in sys.platform: multiprocessing = 'multiprocessing'
elif 'darwin' in sys.platform: multiprocessing = 'multiprocess'
elif 'win' in sys.platform: multiprocessing = 'multiprocessing'
mp = __import__(multiprocessing)

# def give_me_time(secs:float) -> str:

# 	days,seconds = divmod(secs,60*60*24)
# 	hours,seconds = divmod(seconds,60*60)
# 	minutes,seconds = divmod(seconds,60)

# 	if branch_point%10 == 0:
# 		print(f"\t[{branch_point: >{buffer}d}/{total_branch_points}]  Processed 10 branches in {give_me_time(time()-branch_start)}")
# 		print(f"\t\tElapsed Time: {give_me_time(time()-total_start)}")
# 		print(f"\t\tAverage Time Per Branch: {give_me_time(median(times))}")
# 		print(f"\t\tMaximum Branch Time: {give_me_time(max_time)}")
# 		print(f"\t\tMinimum Branch Time: {give_me_time(min_time)}")
# 		print(f"\t\tETA to completion: {give_me_time(median(times)*(total_branch_points-branch_point))}")
# 		branch_start = time()

# 	return f"{int(days)}d {int(hours)}h {int(minutes)}m {int(seconds)}s"

def load_msa(msa_file:str,outdir) -> pd.DataFrame:

	"""
	Read in and store MSA, checking whether length is same for all sequences present
	"""

	print(f"\n\t\tLoading alignment data from {msa_file}")

	data_dir = f"{outdir}/.data"
	pickle_file = f"{data_dir}/MSA.pkl"
	makedirs(data_dir,0o700,exist_ok=True)

	if isfile(pickle_file):
		print(f"\n\t\t\t[N] Found previously processed MSA data")
		print(f"\n\t\t\t[N] Loading previously processed MSA data")
		start = time()
		arrayed_MSA = pd.read_pickle(pickle_file,compression='gzip')
		print(f"\n\t\t\tPreviously processed MSA data loaded in {time()-start:.0f}s")
		return arrayed_MSA

	print(f"\n\t\tReading alignment data from MSA file {msa_file}")

	start = time()
	_,msa = read_MSA(msa_file=msa_file)

	print(f"\n\t\t\tAlignment data loaded in {time()-start:.0f}s")

	msa_keys = msa.keys()

	msa = {x:[i for i in msa[x]] for x in msa_keys}

	len_max = max([len(msa[x]) for x in msa_keys])
	len_min = min([len(msa[x]) for x in msa_keys])

	if len_min != len_max:
		print(f"\n\t[E]  Aligned sequences do not have the same length (min/max: {len_min}/{len_max}).")
		print(f"\n\t[N]  Terminating Two-Entropy calculations.")
		exit()

	pause = time()

	print(f"\n\t\tProcessing alignment data")

	arrayed_MSA = pd.DataFrame({'array':np.array(np.zeros((len(AAs),len_max)) for key in msa_keys)},dtype=object,index=msa_keys)#,columns=[x for x in msa_keys])


	for key in msa_keys:
		for index,res in enumerate(msa[key]):
			arrayed_MSA.loc[key,'array'][AAs.index(res)][index] = 1

	arrayed_MSA.to_pickle(pickle_file,compression='gzip')

	print(f"\n\t\t\tAlignment data prcocessed in {time()-pause:.0f}s")

	return arrayed_MSA

def load_similarity_matrix(matrix:str) -> dict:

	matrix_dir = ("/".join(__file__.split("/")[0:-2]))+"/DB_files"
	labels = []
	sim_mat = {}
	file = None

	if os.path.exists(f"{matrix_dir}/{matrix.upper()}.txt"):
		file = f"{matrix_dir}/{matrix.upper()}.txt"

	elif os.path.exists(matrix):
		file = matrix
	
	if file == None:
		print(f"\n\n\t[E]  {matrix} is not a pre-existing DB, and could not find custom DB file {matrix}\n\n")
		exit()

	with open(file,'r') as IN:
		for line in IN:
			
			line = line.strip()
			
			if line == "":
				continue

			if line[0] == "" or line[0] == "#":
				continue

			if line[0] == "!":
				labels = [x for x in line.split()[1:]]
				sim_mat = {x:{y:0 for y in labels} for x in labels}
				continue

			data = line.split()
			current_label = data[0].strip()

			for index,score in enumerate(data[1:]):
				sim_mat[current_label][labels[index]] = float(score)

		return sim_mat

def process_colummns(msa:pd.DataFrame) -> np.ndarray:

	"""Returns the counts of all unique residues present in all sequences at the i-th index, for all index in the MSA.
	 
	### Input:
	**msa**: _Pandas.DataFrame_
	<br>Multiple-Sequence Alignment where each row is the j-th aligned protein sequence and each column is the i-th residue of the alignment.
	"""
	
	return np.sum(msa.loc[msa.index,'array'],axis=0)

def get_reference_indices(ref:str,msa:pd.DataFrame) -> list:

	if ref not in msa.index:
		print(f"\n\t[W] Reference {ref} is not found in the provided MSA. Skipping to prevent errors\n\n")
		return None
	
	indices = []

	for index,residues in enumerate(msa.loc[ref,'array'].T):
		if residues[-1] == 1:
			continue
		indices.append(index)

	return indices

def calculate_shannon_entropy(msa:pd.DataFrame,weigh:bool=True) -> Tuple[np.ndarray,np.ndarray,np.ndarray]:

	"""
	Calculates the weighted Shannon Entropy for the provided MSA utilizing the proposed equation by Ye, 2008

	#### Input
	**msa**: *pd.DataFrame*
	<br>DataFrame contianing all residue counts for all sequence indexes for all sequences in the MSA</br>

	#### Output
	**E_i**: *np.ndarray*
	<br>Weighted Shannon Entropy for all sequence indexes in the MSA</br>
	**a_i**: *np.ndarray*
	<br>Residue counts for all sequence indexes in the MSA</br>
	**G_i**: *np.ndarray*
	<br>Gap counts for all sequence indexes in the MSA</br>
	"""

	n_j,_ = msa.shape

	_,num_of_MSA_res = msa.loc[msa.index[0],'array'].shape

	## Occurrences of residues at position of i
	residue_counts:dict = process_colummns(msa=msa)

	## Occurrences of AA at position i in the subfamily, minus gaps, normalized by number of proteins in subfamily
	## Array is of size number_of_residues x 20 (number of amino acids)
	a_i:np.array = residue_counts[0:-1]/n_j

	## ln(0) ==> bad things... instead make ln(1) ==> 0 as 0 won't effect the output
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

def generate_sequence_consensus(msa:pd.DataFrame,outdir:str) -> np.ndarray:

	outfile = f"{outdir}/consensus.logo"

	with open(outfile,'w') as OUT:
		OUT.write(f"## MSA_position\tRes:Count;\n")

		residues,msa_length = msa.loc[msa.index[0],'array'].shape

		for msa_index,datum in enumerate(np.sum(msa.loc[msa.index,'array'],axis=0).T):

			OUT.write(f"{msa_index}\t")

			for count,residue in sorted(list(zip(datum,AAs)),key=lambda x: x[0],reverse=True):
				OUT.write(f"{residue}:{int(count)};")
			OUT.write("\n")

	return np.sum(msa.loc[msa.index,'array'],axis=0).T

def write_references(reference_id:str,msa:pd.DataFrame,global_SE:pd.DataFrame,average_SE:np.ndarray,residues:np.ndarray,gaps:np.ndarray,outdir:str,consensus:np.ndarray) -> None:

	print(f"\n\tSummarizing Two Entropy Analysis for reference sequence...")

	if not reference_id: 
		print(f"\n\t\t[N] Reference ID not provided, defaulting to first sequence in MSA [{msa.index[0]}]...")
		reference_id = msa.index[0]

	outfile = f"{outdir}/shannon_entropy.summary"

	reference_seq_to_msa_indexes = get_reference_indices(ref=reference_id,msa=msa)

	with open(outfile,'w') as OUT:
		OUT.write(f"## Sequence_Pos\tMSA_Pos\tResidue\tGlobal_Shannon_Entropy\tAverage_Shannon_Entropy\tSequences_at_pos\tNum_of_gaps\n")

		for res_index,msa_index in enumerate(reference_seq_to_msa_indexes):
			
			global_entropy = global_SE[0][msa_index]
			average_entropy = average_SE[0][msa_index]
			non_gapped_res = int(np.sum(residues[:,msa_index]))
			num_of_gaps = int(gaps[0][msa_index])

			## Residue index
			OUT.write(f"{res_index}")
			## MSA index
			OUT.write(f"\t{msa_index}")
			## Global entropy
			OUT.write(f"\t{global_entropy}")
			## Average entropy
			OUT.write(f"\t{average_entropy}")
			## Residue of reference
			OUT.write(f"\t{AAs[np.where(msa.loc[reference_id,'array'].T[msa_index]==1)[0][0]]}")
			## Consenus
			residue_count = ";".join([f"{residue}:{int(count)}" for count,residue in sorted(list(zip(np.sum(msa.loc[msa.index,'array'],axis=0).T[msa_index],AAs)),key=lambda x: x[0],reverse=True)])
			OUT.write(f"\t{residue_count}")
			## Non-gapped residues
			OUT.write(f"\t{non_gapped_res}")
			## Gapped residues
			OUT.write(f"\t{num_of_gaps}\n")

	return None

def read_subfamilies(subfamilies_file:str) -> dict:

	"""
	Reads subfamily file created by upgma_subfam_grouping.py

	#### Input
	*str* Filepath to subfamilies file


	#### Returns
	*dict* A hash of hashes containing all subfamilies for all branchpoints.
	The first set of keys are the branchpoints and the second is the group numbers for the subfamilies.
	The stored values are all the accessions for the corresponding group at the corresponding
	branchpoint.
	"""

	## Load in subfamily information; headers being keys and values being the subfamily it belongs to
	subfamilies = pd.read_csv(subfamilies_file,header=0,index_col=None,sep=';')

	families = {}

	## Each line represents a branch point; several sequences can be joined at the same branch point
	for index,branch in subfamilies.iterrows():

		families[index] = {}

		## Iterate over all sequences IDs, assigning them to their group numbers
		for key,value in branch.to_dict().items():

			## Initialize subfamliy at the current index
			if value not in families[index].keys():

				families[index][value] = []

			## Add the sequence ID to the proper subfamily
			families[index][value].append(key)

	return families

def calculate_subfamily_shannon_entropy(groupings) -> Tuple[int,np.ndarray]:

	"""Returns the average entropy of all i-th positions within an MSA according to Ye's subfamily weighted averaging.

	### Input:
		**groupings_data**: _tuple_ (branch_index,subfamily_groupings)
		Subfamily groupings, where a key represents the family ID and the value is an array of sequence IDs from the MSA.
	"""

	_,num_of_MSA_res = g_msa.loc[g_msa.index[0],'array'].shape

	se = np.zeros((1,num_of_MSA_res))

	start = time()

	for subfamily_group_number in groupings.keys():

		## Number of proteins in the subfamily
		n_j:int = len(groupings[subfamily_group_number])

		## A subfamily with a single sequence does not provide information on the Two-Entropy approach
		if n_j == 1:
			continue

		E_i,_,_ = calculate_shannon_entropy(msa=g_msa.loc[groupings[subfamily_group_number]])

		se += E_i

	## Return residue-wise entropy normalized by number of subgroups
	return se/len(groupings.keys())

def TEA(subfamilies_file:str,outdir:str,msa:pd.DataFrame,threads:int,mode:str,msa_file:str,file_prefix:str="TEA") -> np.ndarray:

	print("\n\tStarting subfamily entropy calculations...")

	def read_temp_file(temp_file:str) -> dict:

		results = {}

		with open(temp_file,'r') as IN:
			for line in IN:
				line = line.strip()
				if line == "":
					continue
				key,data = line.split("\t")
				results[int(key)] = [float(x) for x in data.split(",")]

		return results

	temp_entropy_file_name = "intermediate_average_entropies.teao"
	if mode == "TEAO":
		subfamilies_file,_ = run_upgma(msa_file=msa_file,outdir=outdir)
		temp_entropy_file_name = "intermediate_average_entropies.teao"
	
	## Default to TEA-O messages
	if mode == "TEA" and not subfamilies_file:
		print("\n\t\t\t[W] Subfamily definition file not provided, defaulting to TEA-O algorithm...")
	if not exists(subfamilies_file):
		print(f"\n\t\t\t[W] Subfamily definition file {subfamilies_file} not found, defaulting to TEA-O algorithm...")
		subfamilies_file,_ = run_upgma(msa_file=msa_file,outdir=outdir)

	subfamilies:dict = read_subfamilies(subfamilies_file)

	## Make the MSA globally accessible
	global g_msa
	g_msa = msa

	temp_entropy_file = f"{outdir}/.data/{temp_entropy_file_name}"

	start = time()

	num_of_seqs,num_of_MSA_res = msa.loc[g_msa.index[0],'array'].shape

	keys = [x for x in subfamilies.keys()]

	results = {}

	if isfile(temp_entropy_file):
		print(f"\n\t\tPrevious subfamily entropies found, resuming at checkpoint...")
		results = read_temp_file(temp_entropy_file)

	queue = mp.Manager().Queue()
	write_temp_process = mp.Process(target=write_temp_results,args=(queue,temp_entropy_file))
	write_temp_process.start()

	with mp.Pool(processes=threads) as executor:
		for subfamily in keys:
			if subfamily in results.keys():
				continue
			executor.apply_async(worker,args=(subfamily,subfamilies[subfamily],queue))#,subfamilies[subfamily],queue,))
		executor.close()
		executor.join()

	queue.put(False)
	write_temp_process.join()

	results = read_temp_file(temp_entropy_file)

	E_i = np.zeros((1,num_of_MSA_res))
	for index in sorted(results.keys()):
		E_i += results[index]
	E_i /= len(results.keys())

	print(f"\n\t\tTEA-O processed {len(keys)} in {time()-start:.0f}s")

	with open(f"{outdir}/{file_prefix}.average_entropy",'w') as OUT:
		for msa_pos,se in enumerate(E_i[0]):
			OUT.write(f"{msa_pos}\t{se:.4f}\n")

	return E_i

def worker(key,subfamily_groupings,queue):
	entropy = calculate_subfamily_shannon_entropy(subfamily_groupings)
	queue.put([key,entropy])
	return entropy

def write_temp_results(queue,temp_entropy_file) -> None:
	with open(temp_entropy_file,'a') as OUT:
		while True:
			data = queue.get()
			if not data:
				break
			key,entropy = data
			entropy = ",".join([f"{x}" for x in entropy[0]])

			OUT.write(f"{key}\t{entropy}\n")
	return None

def run(msa_file:str,subfamilies_file:str,outdir:str,reference_id:str,threads:int,mode:str) -> None:

	print(f"\n\tTwo Entropy Analysis started")
	
	total_start = time()

	global AAs
	AAs = ['A','B','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','X','Z','-']

	# Default settings
	
	makedirs(outdir,mode=0o755,exist_ok=True)

	msa:pd.DataFrame = load_msa(msa_file,outdir=outdir)

	n_seq, n_res = msa.shape

	consensus = generate_sequence_consensus(msa=msa,outdir=outdir)
	
	gE_i,a_i,G_i = superfamily_shannon_entropy(msa=msa,outdir=outdir)

	aE_i = TEA(subfamilies_file=subfamilies_file,outdir=outdir,msa=msa,threads=threads,mode=mode,msa_file=msa_file)

	write_references(reference_id=reference_id,msa=msa,global_SE=gE_i,average_SE=aE_i,residues=a_i,gaps=G_i,outdir=outdir, consensus=consensus)

# def calculate_similarity_score(seq, method='BLOSUM62'):
#     """
#     Using BLOSUM62 matrix. Excluded diagonal and divide by length sequence
#     """
#     aligner = Align.PairwiseAligner()
#     aligner.substitution_matrix = substitution_matrices.load(method)

#     amino_acids = SE.natural_amino_acids()

#     # Remove gaps from sequence
#     # Remove non AA characters from sequence
#     full_length = len(seq)
#     seq = seq.replace('-', '')
#     for char in seq:
#         if char not in amino_acids:
#             seq = seq.replace(char, '')
			
#     num_gaps = full_length - len(seq)
	
#     # Set gap penalty
#     gap_penalty = 0.5 / full_length  # Penalty for each gap. 0.5 is max similarity score

#     # Make matrix
#     M = np.empty((len(seq), len(seq)), dtype='int')

#     # Get upper/lower triangle (including diagonal) and calculate score
#     for i in range(len(seq)):
#         for num, j in enumerate(seq):
#             score = aligner.substitution_matrix[str(seq[i]), str(seq[num])]
#             M[i, num] = float(score)
	
#     # Calculate score by summing values
#     sum_tri = np.sum(M[np.triu_indices(len(seq), k=1)])  # Exclude diagonal
#     sum_diagonal = np.trace(M)  
#     if len(seq) == 1:  # Cannot (and should not) be calculated (Runtimerror)
#         score = np.nan
#     else:
#         score = sum_tri / sum_diagonal  # Divide by diagonal for normalization between residues
#         score = score - (gap_penalty * num_gaps)  # Subtract gap penalty
#         score = round(score / full_length, 2)
	
#     # If score is negative, set to 0
#     if score < 0:
#         score = 0

#     return score


if __name__ == '__main__':

	parser = argparse.ArgumentParser(description='')
	parser.add_argument("-m", "--msa_file", required=True, help="file location of Multiple Sequence Alignment (oneliner format)")
	parser.add_argument("-s", "--subfamily_file", help="file location of subfamily assignment (TXT file)",default=None)
	parser.add_argument("-o", "--outdir", required=True, help="output file location")
	parser.add_argument("-p", "--program_mode", required=False, choices=['TEA', 'TEAO'], default="TEAO", help="run program in TEA or TEAO mode")
	parser.add_argument("-r", "--reference_id", required=False, help="reference ID(s).")
	# parser.add_argument("-i", "--seqidentity", help="Calculate pairwise sequence identity matrix", required=False, default=False)
	parser.add_argument("-x", "--similarity_matrix", help="Calculate scores according to provided matrix", required=False, type=str,default=None)
	# parser.add_argument("-j", "--json", type=str, help="store command line argument as json. Provide folder where to store", required=False)
	parser.add_argument("-t","--threads",type=int,required=False,default=1)
	args = parser.parse_args()

	run(msa_file=args.msa_file,subfamilies_file=args.subfamily_file,outdir=args.outdir,reference_id=args.reference_id,threads=args.threads,mode=args.program_mode)