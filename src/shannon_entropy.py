#!/usr/bin/env python3

import os
import sys
from math import log
from statistics import median
import argparse
import numpy as np
import pandas as pd
from os import makedirs
from os.path import exists
from time import time

multiprocessing = None
if 'linux' in sys.platform: multiprocessing = 'multiprocessing'
elif 'darwin' in sys.platform: multiprocessing = 'multiprocess'
mp = __import__(multiprocessing)

def get_reference_indexs(ref:str,msas:pd.DataFrame) -> list:

	if ref not in msas.index:
		print(f"\n\t[W] Reference {ref} is not found in the provided MSA. Skipping to prevent errors\n\n")
		return 1

	indexs = []
	for index,item in enumerate(msas.loc[ref]):
		if item == '-':
			continue
		indexs.append(index)

	return indexs

def load_msa(msa_file:str) -> pd.DataFrame:

	"""
	Read in and store MSA, checking whether length is same for all sequences present
	"""

	msas = {}

	## Cheap and dirty way to check if all MSA sizes are the same
	size_range = [float('inf'),float('-inf')]
	locus = None
	with open(msa_file,'r') as IN:

		for line in IN:
			
			## Remove trailing spaces and new line characters
			line = line.strip()

			if line == '':
				continue
		
			## Grab Protein_ID
			if line[0] == ">":

				locus = line[1:]

				## Notify of duplicated Protein_ID
				if locus in msas.keys():
				
					print(f"\n\n\t[W]  Locus {locus} has been duplicated in the MSA\n")

				continue

			## Store the MSA string
			msas[locus] = list(line)

			size_range[0] = len(line) if len(line) < size_range[0] else size_range[0]
			size_range[1] = len(line) if len(line) > size_range[1] else size_range[1]

	if size_range[0] != size_range[1]:
		print(f"\n\n\t[E]  Sequences do not have the same length!\n\t\tMin Len: {size_range[0]}\n\t\tMax Len: {size_range[1]}")
		exit()

	msas = pd.DataFrame().from_dict(msas,orient='index')

	return msas

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

def process_colummns(msa:pd.DataFrame) -> dict:

	msas_residues = {}

	_, n_res = msa.shape

	col_data = [msa[x] for x in range(n_res)]

	for res_num,col in enumerate(col_data):

		## Convert column information to array
		residues = col.to_numpy()

		## Get the residue counts
		unique,counts = np.unique(residues,return_counts=True)

		residue_counts = dict(zip(unique,counts))
		msas_residues[res_num] = residue_counts

	return msas_residues

def max_shannon_entropy(num_of_seqs:int) -> float:

	res_prob = (1/num_of_seqs)

	max_sh_entropy = (num_of_seqs * (res_prob * (log(res_prob, 2))))

	return max_sh_entropy

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

def calc_average_entropy(msas:pd.DataFrame) -> list:

	N_j,n_res = msas.shape

	if N_j == 1:
		test = np.zeros(n_res)
		return test

	res_counts = process_colummns(msa=msas)
	AAs = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
	A_j = np.array([[res_counts[res_num][res] if res in res_counts[res_num].keys() else 0 for res in AAs] for res_num in sorted(res_counts.keys())])
	G_j = np.array([res_counts[res_num]['-'] if '-' in res_counts[res_num].keys() else 0 for res_num in sorted(res_counts.keys())])
	W_j = 1/log(N_j) if N_j < 20 else 1/log(20)

	n_j = 1/N_j

	a_j = A_j*n_j
	a_j[a_j == 0] = 1
	
	# print(f"N_j value: {N_j}")
	# print(f"a_j values: {a_j[0]}")
	# print(f"G_j values: {G_j[0]}")
	# print(f"n_j value: {n_j}")
	# print(f"W_j value: {W_j}")

	# print(f"Sum(a_j*ln(a_j): {(-np.sum(((a_j*n_j)*np.log(a_j*n_j))))}")
	# print(f"G_j*n_j*ln(n_j): {(-G_j*n_j*np.log(n_j))[0]}")
	# print(W_j*((-np.sum(((a_j)*np.log(a_j)),axis=1))-(G_j*n_j*np.log(n_j)))[0])

	return W_j*((-np.sum(((a_j)*np.log(a_j)),axis=1))-(G_j*n_j*np.log(n_j)))

def calc_summed_entropy(msas:pd.DataFrame) -> list:
	
	N_j,n_res = msas.shape

	res_counts = process_colummns(msa=msas)
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


def run(args:dict) -> None:

	total_start = time()

	path_msa = args.msapath
	path_subfamilies = args.subfamilypath
	outpath = args.outpath
	reference_ids = args.reference_ids
	threads = args.threads
	mode= args.program_mode
	# sim_mat = args.similarity_matrix

	# Default settings
	shannon_entropy_file = "shannon_entropy.txt"
	tea_shannon_entropy_file = "tea_average_shannon_entropy.txt"
	teao_shannon_entropy_file = "teao_average_shannon_entropy.txt"
	logo_file = "consensus_logo.txt"

	makedirs(outpath,mode=0o755,exist_ok=True)

	# Load and make sure all sequences in the MSA have the same length before starting calculations
	msas = load_msa(path_msa)

	n_seq, n_res = msas.shape
	msa_buffer = len(str(n_res))
	max_sh_entropy = max_shannon_entropy(n_seq)

	start = time()
	results = shannon_entropy(msa=msas)
	min_entropy = min([x[0] for x in results.values()])
	max_entropy = max([x[0] for x in results.values()])

	print(f"Global SE calculated for {n_seq} (nres = {n_res}) in {(time()-start):.2f}s")

	with open(f"{outpath}/{shannon_entropy_file}",'w') as OUT:
		OUT.write(f"## MSA_position\tShannon_Entropy\tFract_Shannon_Entropy\tSequences_at_pos\tNum_of_gaps\n")
		for res_num in sorted(results.keys()):
			sh_entropy,non_gapped_res,num_of_gaps,_ = results[res_num]
			OUT.write(f"{res_num: >{msa_buffer}d}\t{sh_entropy: >.2f}\t{1-((sh_entropy-(min_entropy))/(max_entropy-(min_entropy))): >3.2f}\t{non_gapped_res: >{msa_buffer}}\t{num_of_gaps}\n")

	with open(f"{outpath}/{logo_file}",'w') as OUT:
		OUT.write(f"## MSA_position\tRes:Count;\n")
		for res_num in (sorted(results.keys())):
			_,_,_,residue_counts = results[res_num]
			residues = sorted(residue_counts.keys(),key=lambda x: residue_counts[x],reverse=True)
			OUT.write(f"{res_num: >{msa_buffer}d}\t{residues[0]}:{residue_counts[residues[0]]}")
			for residue in residues[1:]:
				OUT.write(f";{residue}:{residue_counts[residue]}")
			OUT.write("\n")

	for reference in reference_ids:
		indexs = get_reference_indexs(reference,msas)
		res_buffer = len(str(len(indexs)))
		with open(f"{outpath}/{reference}.shannon_entropy.txt",'w') as OUT:
			OUT.write(f"## Sequence_Pos\tMSA_Pos\tResidue\tShannon_Entropy\tFract_Shannon_Entropy\tSequences_at_pos\tNum_of_gaps\n")
			for res_index,msa_index in enumerate(indexs):
				sh_entropy,non_gapped_res,num_of_gaps,_ = results[msa_index]
				OUT.write(f"{res_index: <{res_buffer}d}\t{msa_index: >{msa_buffer}d}\t{msas.loc[reference][msa_index]}\t{sh_entropy: >.2f}\t{(1-(sh_entropy-(min_entropy))/(max_entropy-(min_entropy))): >3.2f}\t{non_gapped_res: >{msa_buffer}}\t{num_of_gaps: >{msa_buffer}}\n")

	if path_subfamilies and exists(path_subfamilies) and (mode == "TEAO"):

		print("Starting subfamily entropy calculations...")

		def give_me_time(secs:float) -> str:

			days,seconds = divmod(secs,60*60*24)
			hours,seconds = divmod(seconds,60*60)
			minutes,seconds = divmod(seconds,60)

			return f"{int(days)}d {int(hours)}h {int(minutes)}m {int(seconds)}s"

		print()

		total_branch_points = 0
		with open(path_subfamilies,'r') as IN:
			total_branch_points = sum(1 for line in IN if line[0] != '#')

		previous_branch_points = 0
		SEs = np.zeros((total_branch_points,n_res))
		temp_file = f"{outpath}/.subfamily_entropies"
		if exists(temp_file):
			with open(temp_file,'r') as TEMP:
				for line_count,line in enumerate(TEMP):
					if line_count < total_branch_points - 1:
						SEs[line_count] = np.fromstring(line,dtype=float,sep=';')
						previous_branch_points += 1

		TEMP = open(f"{outpath}/.subfamily_entropies",'w')

		with open(path_subfamilies,'r') as IN:

			indexes = []

			branch_point = -1
			branch_start = time()
			times = []
			max_time = np.NINF
			min_time = np.PINF
			for line in IN:
				line = line.strip()
				
				if line == "":
					continue
				if line[0] == "#":
					indexes = line.split()[-1].split(';')
					continue

				branch_point += 1
	
				families = {}
				for index,family in enumerate(line.split(";")):
					if family not in families.keys():
						families[family] = []

					families[family].append(indexes[index])

				## Iterate over each subfamily
				start = time()
				m = len(families.keys())

				if branch_point < previous_branch_points:
					TEMP.write(f"{';'.join(map(str, SEs[branch_point]))}\n")
					continue

				results = []
				with mp.Pool(threads) as executor:
					results = executor.map(calc_average_entropy,[msas.loc[family] for family in families.values()])
				results = np.array(results)
				E_i = np.divide(np.sum(results,axis=0),m)

				buffer = len(str(total_branch_points))

				if len(times) == 100:
					times.pop(0)
					times.append(time()-start)
				else:
					times.append(time()-start)

				max_time = times[-1] if times[-1] > max_time else max_time
				min_time = times[-1] if times[-1] < min_time else min_time

				TEMP.write(f"{(';'.join(map(str,E_i)))}\n")

				SEs[branch_point] = E_i
				
				if branch_point%10 == 0:
					print(f"\t[{branch_point: >{buffer}d}/{total_branch_points}]  Processed 10 branches in {give_me_time(time()-branch_start)}")
					print(f"\t\tElapsed Time: {give_me_time(time()-total_start)}")
					print(f"\t\tAverage Time Per Branch: {give_me_time(median(times))}")
					print(f"\t\tMaximum Branch Time: {give_me_time(max_time)}")
					print(f"\t\tMinimum Branch Time: {give_me_time(min_time)}")
					print(f"\t\tETA to completion: {give_me_time(median(times)*(total_branch_points-branch_point))}")
					branch_start = time()

		SEs = np.divide(np.sum(SEs,axis=0),total_branch_points)
		SEs_norm = (SEs - np.min(SEs)) / (np.max(SEs) - np.min(SEs))


		with open(f"{outpath}/{teao_shannon_entropy_file}",'w') as OUT:
			for msa_pos,se in enumerate(SEs_norm):
				OUT.write(f"{msa_pos: <{msa_buffer}}\t{se:.2f}\n")

		print(f"\nTask completed after a total runtime of {give_me_time(time()-total_start)}\n")

	elif path_subfamilies and exists(path_subfamilies) and (mode == "TEA"):

		with open(path_subfamilies,'r') as IN:

					indexes = []

					for line in IN:
						line = line.strip()
						
						if line == "":
							continue
						if line[0] == "#":
							indexes = line.split()[-1].split(';')
							continue
						
						families = {}
						for index,family in enumerate(line.split(";")):
							if family not in families.keys():
								families[family] = []

							families[family].append(indexes[index])
						
						## Iterate over each subfamily
						m = len(families.keys())
						
						results = []
						with mp.Pool(threads) as executor:
							results = executor.map(calc_summed_entropy,[msas.loc[family] for family in families.values()])
						results = np.array(results)
						E_i = np.divide(np.sum(results,axis=0),m)

		with open(f"{outpath}/{tea_shannon_entropy_file}",'w') as OUT:
			for msa_pos,se in enumerate(E_i):
				OUT.write(f"{msa_pos}\t{se:.2f}\n")

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


def main(args):

	run(args)

if __name__ == '__main__':

	parser = argparse.ArgumentParser(description='')
	parser.add_argument("-m", "--msapath", required=True, help="file location of Multiple Sequence Alignment (oneliner format)")
	parser.add_argument("-s", "--subfamilypath", help="file location of subfamily assignment (TXT file)")
	parser.add_argument("-o", "--outpath", required=True, help="output file location")
	parser.add_argument("-p", "--program_mode", required=False, choices=['TEA', 'TEAO'], default="TEAO", help="run program in TEA or TEAO mode")
	parser.add_argument("-r", "--reference_ids", required=False, help="reference ID(s).", nargs="+")
	# parser.add_argument("-i", "--seqidentity", help="Calculate pairwise sequence identity matrix", required=False, default=False)
	parser.add_argument("-x", "--similarity_matrix", help="Calculate scores according to provided matrix", required=False, type=str,default=None)
	# parser.add_argument("-j", "--json", type=str, help="store command line argument as json. Provide folder where to store", required=False)
	parser.add_argument("-t","--threads",type=int,required=False,default=1)
	args = parser.parse_args()

	main(args)