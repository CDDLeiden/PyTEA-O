import os
import time
import textwrap
import xarray as xr
import pandas as pd
import numpy as np
import pathlib
import collections
import typing

from src.utils.multiprocess import Pool
from src.utils.sequence import SequenceUtilities

class MSA:

	_CHAR_TO_INT = {
		'-':0,
		'A':1,
		'B':2,
		'C':3,
		'D':4,
		'E':5,
		'F':6,
		'G':7,
		'H':8,
		'I':9,
		'K':10,
		'L':11,
		'M':12,
		'N':13,
		'P':14,
		'Q':15,
		'R':16,
		'S':17,
		'T':18,
		'V':19,
		'W':20,
		'Y':21,
		'X':22,
		'Z':23
	}

	_INT_TO_CHAR = {
		0:'-',
		1:'A',
		2:'B',
		3:'C',
		4:'D',
		5:'E',
		6:'F',
		7:'G',
		8:'H',
		9:'I',
		10:'K',
		11:'L',
		12:'M',
		13:'N',
		14:'P',
		15:'Q',
		16:'R',
		17:'S',
		18:'T',
		19:'V',
		20:'W',
		21:'Y',
		22:'X',
		23:'Z'
	}

	def __init__(self, msa_file:pathlib.Path,outdir:pathlib.Path,threads:int=1,reference_accession:str|None=None):

		self.msa_file:pathlib.Path = msa_file
		self.outdir:pathlib.Path = outdir
		self.threads:int = threads

		self.__load_msa(msa_file=self.msa_file)

		self.shape = self.msa.shape
		self.sequence_length,self.num_sequences = self.msa.shape
		self.accessions = self.msa.columns
		self.reference_accession = self.__check_reference_accession(reference_accession)
		
		self.descriptors_labels = SequenceUtilities.Sandberg_Zscales['labels']

		self.__residue_counts = None
		self.__distance_matrix = None
		self.__descriptor_deviation = None


	def __load_msa(self,msa_file:pathlib.Path,outdir:pathlib.Path=pathlib.Path('./')) -> None:

		"""
		Reads an MSA from a given file and converts it to a pandas.DataFrame object, where each column is sequence, and each
		row is the i-th position in the MSA

		Inputs:
			MSA_File: File path to MSA, either in FASTA or Clustal format.

		Returns:
			None
		"""

		msa:dict = {}
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

		msa:dict = {x:"".join(msa[x]) for x in msa.keys()}

		msa_lengths:list = [len(x) for x in msa.values()]

		if min(msa_lengths) != max(msa_lengths):
			print(f"\n\t[E]  Aligned sequences do not have the same length (min/max: {min(len(msa_lengths))}/{max(len(msa_lengths))}).")
			print(f"\n\t[N]  Terminating Two-Entropy calculations.")
			exit()

		self.msa = pd.DataFrame({x:[self._CHAR_TO_INT[y] for y in msa[x]] for x in msa.keys()})

	def get_reference_indices(self,ref:str) -> np.array:

		"""
		Returns the indicies of the MSA where the reference has a defined residue

		Parameters:
		ref (str): Accession number of sequence to be used as reference
		msa (Pandas.DataFrame): Multiple-Sequence Alignment where each row is the j-th aligned protein sequence and each column is the i-th residue of the alignment.

		Returns:
		np.array: Indicies where reference has a defined sequence
		"""

		msa = self.msa

		if ref not in msa.columns:
			print(f"\n\t[W] Reference {ref} is not found in the provided MSA. Skipping to prevent errors\n\n")
			return None

		return msa.index[msa[ref]!=0]

	def get_residue_counts(self,accessions:list) -> pd.DataFrame:

		"""
		Returns the counts of all the different residues present in all sequences at the i-th index, for all index in the MSA.
		
		Parameters:
		msa (pandas.DataFrame): Multiple-Sequence Alignment where each row is the j-th aligned protein sequence and each column is the i-th residue of the alignment.

		Returns:
		pandas.DataFrame: Sum of residue counts across all sequences for each column (size = AA[24] x MSA_Length)
		"""

		def __count_residues(row) -> dict:

			counts = collections.Counter([item for item in row if item in self._INT_TO_CHAR])

			return {numeric: counts.get(numeric, 0) for numeric in self._INT_TO_CHAR}
		
		msa:pd.DataFrame = self.msa[accessions]

		## Get residue counts 
		msa_residue_counts = msa.apply(__count_residues,axis=1,result_type='expand')

		msa_df = pd.DataFrame((msa_residue_counts)).fillna(0).astype(int)

		msa_df.columns = [self._INT_TO_CHAR[x] for x in msa_df.columns]

		return msa_df
	
	def get_sequence(self,accession:str|None=None) -> str:

		if accession is None:
			accession = self.reference_accession

		indicies = self.get_reference_indices(accession)

		encoded_sequence = self.msa.loc[indicies][accession]

		return "".join(self.__decode_sequences(encoded_sequence))
	
	def get_concensus_sequence(self) -> str:

		return self.residue_counts.idxmax(axis=1)
	
	def __decode_sequences(self,encoded_residues:list) -> list:

		return [self._INT_TO_CHAR[x] for x in encoded_residues]

	def __calculate_distances(self,sequence_a_loc:str) -> typing.Tuple[str,dict]:

		sequence_a = self.msa[sequence_a_loc]

		distances = {x:0 for x in self.accessions}

		for sequence_b_loc in sorted(distances.keys()):
			if sequence_b_loc == sequence_a_loc:
				distances[sequence_b_loc] = 0
				return sequence_a_loc,distances
			distances[sequence_b_loc] = np.sum(sequence_a!=self.msa[sequence_b_loc])/2

		return sequence_a_loc,distances

	def __calculate_sequence_difference(self) -> pd.DataFrame:

		keys:list = sorted(self.accessions)

		distance_matrix:pd.DataFrame = pd.DataFrame(np.zeros([len(keys),len(keys)]),columns=keys,index=keys)

		results:list
		with Pool(self.threads) as executor:
			results = executor.map(self.__calculate_distances,keys)

		for sequence_loc_a,result in results:
			distance_matrix.loc[sequence_loc_a] = result

		distance_matrix += distance_matrix.T

		return distance_matrix
	
	def __apply_sequence_descriptors(self):

		## WxR matrix
		descriptors:pd.DataFrame = pd.DataFrame.from_dict(SequenceUtilities.Sandberg_Zscales['values'])

		descriptors.index = index=SequenceUtilities.Sandberg_Zscales['labels']

		## PxR
		residue_counts:pd.DataFrame = self.residue_counts[self.residue_counts.columns.intersection(descriptors.columns)]

		## 1xRxW
		weights = descriptors.values.T[np.newaxis:,:]
		## PxRx1
		counts = residue_counts.values[:,:,np.newaxis]

		## Px1
		total_res = residue_counts.sum(axis=1).values[:,np.newaxis]

		## Px3
		weighted_mean = np.sum(counts*weights,axis=1)/total_res

		diff_squared = (weights[np.newaxis,:,:]-weighted_mean[:,np.newaxis,:])**2
		weighted_variance = np.sum(counts*diff_squared,axis=1)/total_res

		standard_dev = np.sqrt(weighted_variance)

		max_dev = standard_dev.max().max()

		normalized_dev = pd.DataFrame(standard_dev/max_dev,columns=descriptors.index)

		return normalized_dev

	def __check_reference_accession(self,accession=str|None) -> str:
		reference_accession:str
		if accession in self.accessions:
			reference_accession = accession
		else:
			reference_accession = self.accessions[0]
		return reference_accession
	
	@property
	def distance_matrix(self) -> pd.DataFrame:
		if self.__distance_matrix is None:
			self.__distance_matrix = self.__calculate_sequence_difference()
		return self.__distance_matrix
	
	@property
	def residue_counts(self) -> pd.DataFrame:
		if self.__residue_counts is None:
			self.__residue_counts = self.get_residue_counts(accessions=self.msa.columns)
		return self.__residue_counts

	@property
	def descriptor_deviation(self) -> pd.DataFrame:
		if self.__descriptor_deviation is None:
			self.__descriptor_deviation = self.__apply_sequence_descriptors()
		return self.__descriptor_deviation


