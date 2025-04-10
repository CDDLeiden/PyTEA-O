#!/usr/bin/env python3

import os
import sys
import argparse
import numpy as np
import pandas as pd
from time import time


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


def read_MSA(MSA_File:str=None) -> dict:

	"""Reads an MSA from a given file, and converts it to a dictionary object, with loci and sequences as keys and
	values, repectively.

	Inputs:
		MSA_File: File path to MSA, either in FASTA or Clustal format.

	Returns:
		MSA: Dictionary object containig MSAs
	"""

	start = time()

	msa = {}
	with open (MSA_File, 'r') as f:

		loci = None

		for i, line in enumerate(f):

			line = line.strip()
			if line == '':
				continue
	
			split_line = line.split()

			if split_line[0][0] == ">":
				loci = split_line[0][1:]
				msa[loci] = []
				continue
			elif len(split_line) == 2:
				loci = split_line[0]
				line = split_line[1]

			msa[loci].append(line)

	msa = {x:set([(i,char) for i,char in enumerate("".join(msa[x]))]) for x in list(msa.keys())}

	return msa

def generate_distance_matrix(msa:dict) -> pd.DataFrame:

	"""Returns an NxN numpy matrix with values set to -1.

	Inputs:
		msa [dict]: Dictionary of MSAs, where keys are the loci and values are the sequences

	Returns:
		np.array: NxN array with each element being set to -1
	"""

	keys:dict.keys = sorted(msa.keys())

	return pd.DataFrame(np.full([len(msa),len(msa)],np.inf),columns=keys,index=keys)

def distance(sequence_a_loc):

	sequence_a = sequences[sequence_a_loc]

	distances = {}

	for sequence_b_loc in sorted(sequences.keys()):
		if sequence_b_loc == sequence_a_loc:
			break
		distances[sequence_b_loc] = len(sequence_a.difference(sequences[sequence_b_loc]))

	return [sequence_a_loc,distances]

def calculate_sequence_difference(msa:dict) -> pd.DataFrame:

	"""Calculates the distances between sequences in an MSA

	Returns:
		pandas.Dataframe: NxN lower triangular matrix
	"""

	start = time()

	distance_matrix = generate_distance_matrix(msa)

	keys:dict.keys = sorted(msa.keys())
	global sequences
	sequences = msa

	with mp.Pool(mp.cpu_count()) as executor:
		results = executor.map(distance,keys)

	print(f"Calculated {len(keys)**2/2} distances in {int(time()-start)}s")


	exit()
	return distance_matrix

def grouping(msa:dict, outfile:str):

	distance_matrix:pd.DataFrame = calculate_sequence_difference(msa)
	
	keys = list(msa.keys())

	original_distance_matrix = distance_matrix.copy()
	# location = []
	# for i in range(0,len(distance_matrix)):
	# 	location.append([i])

	# print(location)

	# with open(outfile, "w") as file:
	# 	first_line = "## " + ";".join(keys) + "\n"
	# 	file.write(first_line)

	tree = {
		0 : keys
	}

	# start_time = datetime.now()
	while(len(distance_matrix)>2):
		# row = 0
		# m = len(seqMatrix)
		flattend_index =  np.argmin(distance_matrix)

		row = keys[flattend_index//distance_matrix.shape[1]]
		col = keys[flattend_index%distance_matrix.shape[1]]

		print(row,col,distance_matrix.loc[row,col])

		exit()
		# while(index-(m-row-1)>=0): 
		# 	index -= (m-row)-1
		# 	row += 1
		# column = 1
		# while(index-1>=0):
		# 	index -= 1
		# 	column += 1
		# column += row
	
		#make newick format
		key1 = keys[row] # Name in the (row)th value of key list
		key2 = keys[column] # Name in the (column)th value of the key list
		length = seqMatrix[row, column]/2
		
		new_key = "{},{}".format(key1, key2)

		del keys[column]
		del keys[row]

		keys.insert(0,new_key)
		'''
		this is for teao
		'''
		teao_list = keys.copy()
		for i in range(len(teao_list)):
			if ',' in teao_list[i]:
				teao_list[i] = teao_list[i].split(',')
		tree[length] = teao_list

		# Make new array and update it
		arr2 = seqMatrix
		arr2 = np.delete(arr2, column, axis = 0)
		arr2 = np.delete(arr2, column, axis = 1)

		arr2 = np.delete(arr2, row, axis = 0)
		arr2 = np.delete(arr2, row, axis = 1)
		
		arr2 = np.insert(arr2, 0, -1, axis = 1) # (arr, obj, values, axis)
		arr2 = np.insert(arr2, 0, -1, axis = 0)


		new_name = location[row]+location[column]
		location.insert(0,new_name)

		del location[column+1]
		del location[row+1]

		for l in range(1, len(arr2)):
			sum1 = 0
			for m in range(0, len(location[0])): # [abc]
				for k in range(0, len(location[l])): # [de],[f],[g]
					if (location[0])[m] > (location[l])[k]:
						sum1 += original_seqM[(location[l])[k], (location[0])[m]]
					else:
						sum1 += original_seqM[(location[0])[m],(location[l])[k]]
					arr2[0,l] = sum1/((len(location[0]))*(len(location[l])))

		seqMatrix = arr2
	end_time = datetime.now()
	print(end_time-start_time)
	
	# Final newick
	row = 0
	column = 1
	tree[seqMatrix[row,column]/2] = [tree[0]]

	# level = 0
	tree_red = {}
	tree_content = {}

	for j, value in enumerate(tree.values()):
		for i in range(len(value)):
			if type(value[i])==list:
				# print(value[i])
				for k in range(len(value[i])):
					# print(value[i][k])
					if type(value[i][k])==str:
						tree_content[value[i][k]] = hex(i)
					elif type(value[i][k])==list:  # Duplicate sequences
						print(value[i][k])
						sys.exit()
						tree_content[str(value[i][k])] = hex(i)
			else:
				# print(value[i])
				tree_content[value[i]] = hex(i)
		content = tree_content.copy()
		tree_red[j] = content
		grouping_list = list(content.values())
		grouping_line = ";".join(grouping_list)
		with open(outfile, 'a') as file:
			file.write(grouping_line+"\n")

def main(args):

	start = time()

	msa = read_MSA(MSA_File=args.msa)

	# Create output filename
	basename = os.path.splitext(args.msa)[0]
	outname = os.path.join(basename + "_grouping.txt")
	grouping(msa,outname)
	exit()



if __name__ == '__main__':

	parser = argparse.ArgumentParser(description="""Grouping subfamilies based on UPGMA clustering. Input is a Multiple Sequence Alignment (CLUSTAL format)
													and the output is a text file containing the subfamily grouping. The output file is named as 
													<input_filename>_grouping.txt""")
	parser.add_argument("-m", "--msa", required=True, help="file location of Multiple Sequence Alignment (CLUSTAL format)")
	args = parser.parse_args()

	main(args)