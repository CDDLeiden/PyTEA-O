#!/usr/bin/env python3

from argparse import ArgumentParser
import os
from textwrap import wrap
import pandas as pd
import numpy as np
from src.utils.sequence import SequenceUtilities
from src.utils.general import create_directory
import time

def get_consensus_sequence(msa:pd.DataFrame,outdir:str) -> np.ndarray:

	outfile = f"{outdir}/consensus.logo"
	consensus = []

	with open(outfile,'w') as OUT:
		
		OUT.write(f"## MSA_position\tRes:Count;\n")

		for msa_index,datum in enumerate(np.sum(msa.loc[msa.index,'array'],axis=0).T):

			OUT.write(f"{msa_index}\t")

			for count,residue in sorted(list(zip(datum,SequenceUtilities.AAs)),key=lambda x: x[0],reverse=True):
				consensus.append(f"{residue}:{int(count)}")
				OUT.write(f"{consensus[-1]};")
			OUT.write("\n")

	return ";".join(consensus)

def read_MSA(msa_file:str=None) -> tuple[dict,str]:

	"""Reads an MSA from a given file, and converts it to a dictionary object, with loci and sequences as keys and
	values, repectively.

	Inputs:
		MSA_File: File path to MSA, either in FASTA or Clustal format.

	Returns:
		MSA: Dictionary object containig MSAs
	"""

	file_prefix:str = msa_file.split("/")[-1].split(".")[0]

	msa = {}
	with open (msa_file, 'r') as MSA_FILE:

		loci = None

		for i, line in enumerate(MSA_FILE):

			line = line.strip()
			if line == '':
				continue
	
			split_line = line.split()

			if split_line[0][0] == ">":
				loci = split_line[0][1:].split(".")[0]
				msa[loci] = []
				continue
			elif len(split_line) == 2:
				loci = split_line[0]
				line = split_line[1]

			msa[loci].append(line.upper())

	msa = {x:"".join(msa[x]) for x in list(msa.keys())}

	return file_prefix,msa

def to_one_liner(aligned_sequences:dict,output:str,reference:str) -> None:

	with open(output,'w') as OUT:

		if reference:

			OUT.write(f">{reference}\n{aligned_sequences[reference]}\n")

		for sequence in sorted(aligned_sequences.keys()):

			if sequence == "reference":
				next

			OUT.write(f">{sequence}\n{aligned_sequences[sequence]}\n")

	return None

def extract(aligned_sequences:dict,outdir=str) -> None:

	fasta_dir = f"{outdir}/FASTAs" 
	create_directory(f"{outdir}/FASTAs")

	for locus in aligned_sequences.keys():

		with open(f"{fasta_dir}/{locus}.fasta",'w') as FASTA:

			sequence = "\n".join(wrap(aligned_sequences[locus].replace("-",""),60))

			FASTA.write(f">{locus}\n{sequence}\n")

	return None

def convert(aligned_sequences:dict,align_type:str,outdir:str,file_prefix:str,reference:str) -> None:

	if reference and reference not in aligned_sequences.keys():

		print(f"\n[W]  Reference locus - {reference} - was not found in the alignment set. Defaulting to normal functionality.\n")

		reference = None

	if align_type == "one_liner":
		
		to_one_liner(aligned_sequences=aligned_sequences,output=f"{outdir}/{file_prefix}.ola",reference=reference)
	
	if align_type == "clustal":
		...
	
	if align_type == "fasta":
		...

	return

def load_msa(msa_file:str,outdir) -> pd.DataFrame:

	"""
	Read in and store MSA, checking whether length is same for all sequences present
	"""

	print(f"\n\t\tLoading alignment data from {msa_file}")

	AAs = SequenceUtilities.AAs

	data_dir = f"{outdir}/.data"
	pickle_file = f"{data_dir}/MSA.pkl"
	create_directory(data_dir)

	if os.path.isfile(pickle_file):
		print(f"\n\t\t\t[N] Found previously processed MSA data")
		print(f"\n\t\t\t[N] Loading previously processed MSA data")
		start = time.time()
		arrayed_MSA = pd.read_pickle(pickle_file,compression='gzip')
		print(f"\n\t\t\tPreviously processed MSA data loaded in {time.time()-start:.0f}s")
		return arrayed_MSA

	print(f"\n\t\tReading alignment data from MSA file {msa_file}")

	start = time.time()
	msa = read_MSA(msa_file=msa_file)

	print(f"\n\t\t\tAlignment data loaded in {time.time()-start:.0f}s")

	msa_keys = [x for x in msa.keys()]

	msa = {x:[i for i in msa[x]] for x in msa_keys}

	len_max = max([len(msa[x]) for x in msa_keys])
	len_min = min([len(msa[x]) for x in msa_keys])

	if len_min != len_max:
		print(f"\n\t[E]  Aligned sequences do not have the same length (min/max: {len_min}/{len_max}).")
		print(f"\n\t[N]  Terminating Two-Entropy calculations.")
		exit()

	pause = time.time()

	print(f"\n\t\tProcessing alignment data")

	## Create a database where each sequence has a flattened matrix that is the size of the alignment length x the number of recognized amino acids
	arrayed_MSA = pd.DataFrame({'array':np.array(np.zeros((len(AAs),len_max)) for key in msa_keys)},dtype=object,index=msa_keys)

	## Iterate over every sequence
	for key in msa_keys:

		## Iterate through the entire sequence
		for index,res in enumerate(msa[key]):
			
			## Assign the residue at the current position a value of 1
			arrayed_MSA.loc[key,'array'][AAs.index(res)][index] = 1

	# arrayed_MSA.to_pickle(pickle_file,compression='gzip')

	print(f"\n\t\t\tAlignment data prcocessed in {time.time()-pause:.0f}s")

	return arrayed_MSA

def read_MSA(msa_file:str=None) -> dict:

	"""Reads an MSA from a given file, and converts it to a dictionary object, with loci and sequences as keys and
	values, repectively.

	Inputs:
		MSA_File: File path to MSA, either in FASTA or Clustal format.

	Returns:
		MSA: Dictionary object containig MSAs
	"""

	msa = {}
	with open (msa_file, 'r') as MSA_FILE:

		loci = None

		for i, line in enumerate(MSA_FILE):

			line = line.strip()
			if line == '':
				continue
	
			split_line = line.split()

			if split_line[0][0] == ">":
				loci = split_line[0][1:].split(".")[0]
				msa[loci] = []
				continue
			elif len(split_line) == 2:
				loci = split_line[0]
				line = split_line[1]

			msa[loci].append(line.upper())

	return {x:"".join(msa[x]) for x in list(msa.keys())}

def get_reference_indices(ref:str,msa:pd.DataFrame) -> np.array:

	"""
	Returns the indicies of the MSA where the reference has a defined residue

	Parameters:
	ref (str): Accession number of sequence to be used as reference
	msa (Pandas.DataFrame): Multiple-Sequence Alignment where each row is the j-th aligned protein sequence and each column is the i-th residue of the alignment.

	Returns:
	np.array: Indicies where reference has a defined sequence
	"""

	if ref not in msa.index:
		print(f"\n\t[W] Reference {ref} is not found in the provided MSA. Skipping to prevent errors\n\n")
		return None

	return np.where((msa.loc[ref,'array'].T != 1)[:,-1])[0]

def get_residue_counts(msa:pd.DataFrame) -> np.ndarray:

	"""
	Returns the counts of all the different residues present in all sequences at the i-th index, for all index in the MSA.
	 
	Parameters:
	msa (Pandas.DataFrame): Multiple-Sequence Alignment where each row is the j-th aligned protein sequence and each column is the i-th residue of the alignment.

	Returns:
	np.ndarray: Sum of residue counts across all sequences for each column (size = AA[24] x MSA_Length)
	"""
	
	return np.sum(msa['array'],axis=0)

def run(msa_file:str,outdir:str,extract_fastas:str,desired_alignment_type:str,reference:str) -> int:

	create_directory(outdir)

	file_prefix:str
	alignment_data:dict

	file_prefix,alignment_data = read_MSA(msa_file=msa_file)

	convert( aligned_sequences=alignment_data,
			 align_type=desired_alignment_type,
			 outdir=outdir,
			 file_prefix=file_prefix,
			 reference=reference
	)

	if extract_fastas: extract(aligned_sequences=alignment_data,outdir=outdir)

	return 0

if __name__ == "__main__":

	parser = ArgumentParser()

	parser.add_argument("-m","--msa_file",required=True,type=str)
	parser.add_argument("-o","--outdir",default="alignment_conversion")
	parser.add_argument("-x","--extract_fastas",default=False,action='store_true')
	parser.add_argument("-a","--desired_alignment_type",default=["one_liner"],choices=["one_liner,clustal,fasta"])
	parser.add_argument("-r","--reference",default=None,type=str)

	run(parser.parse_known_args()[0])