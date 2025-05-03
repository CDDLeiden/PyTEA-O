from os import makedirs,path
from pandas import DataFrame,read_pickle,to_pickle
from time import time
from numpy import zeros,array

def load_msa(msa_file:str,outdir) -> DataFrame:

	"""
	Read in and store MSA, checking whether length is same for all sequences present
	"""

	print(f"\n\t\tLoading alignment data from {msa_file}")

	AAs = ['A','B','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','X','Z','-']

	data_dir = f"{outdir}/.data"
	pickle_file = f"{data_dir}/MSA.pkl"
	makedirs(data_dir,0o700,exist_ok=True)

	if path.isfile(pickle_file):
		print(f"\n\t\t\t[N] Found previously processed MSA data")
		print(f"\n\t\t\t[N] Loading previously processed MSA data")
		start = time()
		arrayed_MSA = read_pickle(pickle_file,compression='gzip')
		print(f"\n\t\t\tPreviously processed MSA data loaded in {time()-start:.0f}s")
		return arrayed_MSA

	print(f"\n\t\tReading alignment data from MSA file {msa_file}")

	start = time()
	msa = read_MSA(msa_file=msa_file)

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

	arrayed_MSA = DataFrame({'array':array(zeros((len(AAs),len_max)) for key in msa_keys)},dtype=object,index=msa_keys)#,columns=[x for x in msa_keys])

	for key in msa_keys:
		for index,res in enumerate(msa[key]):
			arrayed_MSA.loc[key,'array'][AAs.index(res)][index] = 1

	arrayed_MSA.to_pickle(pickle_file,compression='gzip')

	print(f"\n\t\t\tAlignment data prcocessed in {time()-pause:.0f}s")

	return arrayed_MSA

def read_MSA(msa_file:str=None) -> tuple[dict,str]:

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

	msa = {x:"".join(msa[x]) for x in list(msa.keys())}

	return msa