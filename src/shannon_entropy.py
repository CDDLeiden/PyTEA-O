#!/usr/bin/env python3

import os
import sys
from math import log
from statistics import median
import argparse
import numpy as np
from numpy.typing import NDArray
import pandas as pd
from os import makedirs
from os.path import exists
from time import time
from convert_alignment import read_MSA
from typing import Tuple

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

global AAs
AAs = ['A','B','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','X','Z','-']

def load_msa(msa_file:str) -> pd.DataFrame:

	"""
	Read in and store MSA, checking whether length is same for all sequences present
	"""

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

def calc_shannon_entropy(res_count:dict) -> float:

	"""
	Calculate Shannon's Entropy per column of an alignment (H=-\\sum_{i=1}^{M} P_i\\,log_2\\,P_i)
	"""

	num_of_seqs = sum([x for x in res_count.values()])

	entropy = 0
	
	# Number of residues in column
	for residue in res_count.keys():
		
		coeff = 1
		num_of_res = res_count[residue]

		if residue == '-':
			coeff = res_count[residue]
			num_of_res = 1

		prob = num_of_res / num_of_seqs

		entropy += (coeff * (prob * (log(prob, 2))))

	return entropy

def calc_summed_entropy(branchpoint:dict) -> list:
	
	N_j,n_res = g_msa.shape

	res_counts = process_colummns(msa=msa)
	AAs = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
	F_ia = {}
	E_i = []
	for position, counts in res_counts.items():
		F_ia[position] = []
		entropy = 0
		for aa, count in counts.items():
			if aa == '-':
				for i in range(count):
					F_ia[position].append(1/N_j)
			elif aa in AAs:
				F_ia[position].append(count/N_j)
		entropy += -1 * np.sum(F_ia[position]*np.log(F_ia[position]))
		E_i.append(entropy)
	return E_i

def get_reference_indexs(ref:str,msa:pd.DataFrame) -> list:

	if ref not in msa.index:
		print(f"\n\t[W] Reference {ref} is not found in the provided MSA. Skipping to prevent errors\n\n")
		return None

	indexes = []
	for index,item in enumerate(msa.loc[ref]):
		if item == '-':
			continue
		indexes.append(index)

	return indexes

def shannon_entropy(msa:pd.DataFrame) -> dict:

	n_seq, n_res = msa.shape

	res_counts = process_colummns(msa)

	results = {}
	for res_num in res_counts.keys():
		entropy = calc_shannon_entropy(res_count=res_counts[res_num])
		num_of_gaps = 0 if '-' not in res_counts[res_num].keys() else res_counts[res_num]['-']
		num_of_non_gaps = n_seq - num_of_gaps
		results[res_num] = [entropy,num_of_non_gaps,num_of_gaps,res_counts[res_num]]

	return results

def write_shannon_entropy(results:dict,msa_buffer:int,outfile:str) -> None:

	min_entropy = min([x[0] for x in results.values()])
	max_entropy = max([x[0] for x in results.values()])

	with open(outfile,'w') as OUT:
		OUT.write(f"## MSA_position\tShannon_Entropy\tFract_Shannon_Entropy\tSequences_at_pos\tNum_of_gaps\n")
		for res_num in sorted(results.keys()):
			
			sh_entropy,non_gapped_res,num_of_gaps,_ = results[res_num]
			
			## MSA Position
			OUT.write(f"{res_num: >{msa_buffer}d}")

			## Shannon Entropy
			OUT.write(f"\t{sh_entropy: >.2f}")

			## Fraction Shannon Entropy
			OUT.write(f"\t{1-((sh_entropy-(min_entropy))/(max_entropy-(min_entropy))): >3.2f}")

			## Sequences at position
			OUT.write(f"\t{non_gapped_res: >{msa_buffer}}")

			## Number of gaps
			OUT.write(f"\t{num_of_gaps}\n")

	return None

def write_logo(results:dict,msa_buffer:int,outfile:str) -> None:

	with open(outfile,'w') as OUT:
		OUT.write(f"## MSA_position\tRes:Count;\n")
		for res_num in (sorted(results.keys())):
			_,_,_,residue_counts = results[res_num]
			residues = sorted(residue_counts.keys(),key=lambda x: residue_counts[x],reverse=True)
			OUT.write(f"{res_num: >{msa_buffer}d}\t{residues[0]}:{residue_counts[residues[0]]}")
			for residue in residues[1:]:
				OUT.write(f";{residue}:{residue_counts[residue]}")
			OUT.write("\n")

	return None

def write_references(reference_id:str,msa:pd.DataFrame,results:dict,outdir:str,msa_buffer:int) -> None:

	min_entropy = min([x[0] for x in results.values()])
	max_entropy = max([x[0] for x in results.values()])

	res_buffer = len(str(msa.shape[1]))
	outfile = f"{outdir}/{reference_id}.shannon_entropy.txt"

	reference_seq_to_msa_indexes = get_reference_indexs(ref=reference_id,msa=msa)

	with open(outfile,'w') as OUT:
		OUT.write(f"## Sequence_Pos\tMSA_Pos\tResidue\tShannon_Entropy\tFract_Shannon_Entropy\tSequences_at_pos\tNum_of_gaps\n")
		res_index = 0
		for res_index,msa_index in enumerate(reference_seq_to_msa_indexes):
			sh_entropy,non_gapped_res,num_of_gaps,_ = results[msa_index]
			OUT.write(f"{res_index: <{res_buffer}d}")
			OUT.write(f"\t{msa_index: >{msa_buffer}d}")
			OUT.write(f"\t{msa.loc[reference_id][msa_index]}")
			OUT.write(f"\t{sh_entropy: >.2f}")
			OUT.write(f"\t{(1-(sh_entropy-(min_entropy))/(max_entropy-(min_entropy))): >3.2f}")
			OUT.write(f"\t{non_gapped_res: >{msa_buffer}}")
			OUT.write(f"\t{num_of_gaps: >{msa_buffer}}\n")
			res_index += 1

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

def calc_average_entropy(index,groupings) -> Tuple[int,NDArray]:

	"""Returns the average entropy of all i-th positions within an MSA according to Ye's subfamily weighted averaging.

	### Input:
		**groupings_data**: _tuple_ (branch_index,subfamily_groupings)
		Subfamily groupings, where a key represents the family ID and the value is an array of sequence IDs from the MSA.
	"""

	_,num_of_MSA_res = g_msa.loc[g_msa.index[0],'array'].shape

	E_i = np.zeros((1,num_of_MSA_res))

	start = time()

	print(f"\n\t\tStarting calculations for Branch {index}. Averaging across {len(groupings.keys())} subfamilies")

	for subfamily_group_number in groupings.keys():

		## Number of proteins in the subfamily
		n_j:int = len(groupings[subfamily_group_number])

		## A subfamily with a single sequence does not provide information on the Two-Entropy approach
		if n_j == 1:
			continue

		## Occurrences of residues at position of i
		res_i:dict = process_colummns(msa=g_msa.loc[groupings[subfamily_group_number]])

		## Occurrences of AA at position i in the subfamily, minus gaps, normalized by number of proteins in subfamily
		## Array is of size number_of_residues x 20 (number of amino acids)
		a_i:np.array = res_i[0:-1]/n_j

		## ln(0) ==> bad things... instead make ln(1) ==> 0 as 0 won't effect the output
		a_i[a_i==0] = 1

		## Number of gaps at position in the subfamily
		G_i:np.array = res_i[-1]

		G_i = G_i.reshape(1,num_of_MSA_res)

		## Normalization factor
		W_j:float = 1/np.log(n_j) if n_j < 20 else 1/np.log(20)

		## Add weighted subfamily entropy to total entropy
		E_i += W_j*((-np.sum(a_i*np.log(a_i),axis=0)).reshape(1,num_of_MSA_res)-G_i*(1/n_j)*np.log(1/n_j))

	print(f"\n\t\tFinished processing Branch {index} in {time()-start:.0f}s")

	## Return residue-wise entropy normalized by number of subgroups
	return index,E_i/len(groupings.keys())

def TEAO(subfamilies_file:str,outdir:str,msa:pd.DataFrame,threads:int) -> None:

	## Make the MSA globally accessible
	global g_msa
	g_msa = msa

	if not exists(subfamilies_file):

		return None

	print("\n\tStarting subfamily entropy calculations...")

	num_of_seqs,num_of_MSA_res = msa.loc[g_msa.index[0],'array'].shape

	msa_buffer = len(str(num_of_seqs))

	subfamilies:dict = read_subfamilies(subfamilies_file)

	keys = [x for x in subfamilies.keys()]

	results = []
	
	for key in keys:
		results.append(calc_average_entropy(key,subfamilies[key]))

	buffer = len(str(len(keys)))

	se = np.zeros((1,num_of_MSA_res))
	for index,result in sorted(results,key=lambda x: x[0]):
		se += result
	se /= len(results)

	with open(f"{outdir}/teao_average_shannon_entropy.txt",'w') as OUT:
		for msa_pos,se in enumerate(se[0]):
			OUT.write(f"{msa_pos: <{msa_buffer}}\t{se:.4f}\n")

	return None

def run(path_msa:str,subfamilies_file:str,outdir:str,reference_id:str,threads:int,mode:str) -> None:

	total_start = time()

	print(f"\n\tTwo Entropy Analysis started")

	# Default settings
	shannon_entropy_file = "shannon_entropy.txt"
	tea_shannon_entropy_file = "tea_average_shannon_entropy.txt"
	logo_file = "consensus_logo.txt"

	makedirs(outdir,mode=0o755,exist_ok=True)

	msa = load_msa(path_msa)

	n_seq, n_res = msa.shape

	msa_buffer = len(str(n_res))

	start = time()
	# results = shannon_entropy(msa=msa)

	# print(f"Global SE calculated for {n_seq} (nres = {n_res}) in {(time()-start):.2f}s")

	# write_shannon_entropy(results=results,msa_buffer=msa_buffer,outfile=f"{outdir}/{shannon_entropy_file}")
	
	# write_logo(results=results,msa_buffer=msa_buffer,outfile=f"{outdir}/{logo_file}")

	# write_references(reference_id=reference_id,msa=msa,results=results,outdir=outdir,msa_buffer=msa_buffer)


	print("\n\tRunning TEA-O...")
	TEAO(subfamilies_file=subfamilies_file,outdir=outdir,msa=msa,threads=threads)

	# if subfamilies_file and exists(subfamilies_file) and (mode == "TEA"):

	# 	with open(subfamilies_file,'r') as IN:

	# 				indexes = []

	# 				for line in IN:
	# 					line = line.strip()
						
	# 					if line == "":
	# 						continue
	# 					if line[0] == "#":
	# 						indexes = line.split()[-1].split(';')
	# 						continue
						
	# 					families = {}
	# 					for index,family in enumerate(line.split(";")):
	# 						if family not in families.keys():
	# 							families[family] = []

	# 						families[family].append(indexes[index])
						
	# 					## Iterate over each subfamily
	# 					m = len(families.keys())
						
	# 					results = []
	# 					with mp.Pool(threads) as executor:
	# 						results = executor.map(calc_summed_entropy,[msa.loc[family] for family in families.values()])
	# 					results = np.array(results)
	# 					E_i = np.divide(np.sum(results,axis=0),m)
	# 					E_i_norm = (E_i - np.min(E_i)) / (np.max(E_i) - np.min(E_i))


		# with open(f"{outdir}/{tea_shannon_entropy_file}",'w') as OUT:
		# 	for msa_pos,se in enumerate(E_i_norm):
		# 		OUT.write(f"{msa_pos}\t{se:.2f}\n")

	# load_similarity_matrix(matrix=sim_mat)

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
	parser.add_argument("-m", "--msapath", required=True, help="file location of Multiple Sequence Alignment (oneliner format)")
	parser.add_argument("-s", "--subfamilypath", help="file location of subfamily assignment (TXT file)")
	parser.add_argument("-o", "--outdir", required=True, help="output file location")
	parser.add_argument("-p", "--program_mode", required=False, choices=['TEA', 'TEAO'], default="TEAO", help="run program in TEA or TEAO mode")
	parser.add_argument("-r", "--reference_id", required=False, help="reference ID(s).")
	# parser.add_argument("-i", "--seqidentity", help="Calculate pairwise sequence identity matrix", required=False, default=False)
	parser.add_argument("-x", "--similarity_matrix", help="Calculate scores according to provided matrix", required=False, type=str,default=None)
	# parser.add_argument("-j", "--json", type=str, help="store command line argument as json. Provide folder where to store", required=False)
	parser.add_argument("-t","--threads",type=int,required=False,default=1)
	args = parser.parse_args()

	run(path_msa=args.msapath,subfamilies_file=args.subfamilypath,outdir=args.outdir,reference_id=args.reference_id,threads=args.threads,mode=args.program_mode)