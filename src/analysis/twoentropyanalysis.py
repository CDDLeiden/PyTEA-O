#!/usr/bin/env python3

import numpy as np
import pandas as pd
from time import time
import typing
import heapq
import pathlib

from src.utils.msa import MSA
from src.analysis.phylotree import PhyloTree
from src.utils.general import valid_directory

class TwoEntropyAnalysis:

	def __init__(self,msa:MSA,tree:PhyloTree,threads:int=1,outdir:str="./TEA"):

		self.msa = msa
		self.tree = tree
		self.threads = threads
		self.outdir = valid_directory(outdir)
		self.__average_entropy = None
		self.__global_entropy = None

	def __calculate_shannon_entropy(self,node_id:int,weigh:bool=True):

		"""
		Calculates the weighted Shannon Entropy for the provided MSA utilizing the proposed equation by Ye, 2008

		Parameters
		msa (pd.DataFrame): DataFrame containing all residue counts for all sequence indexes for all sequences in the MSA

		Output
		E_i (np.ndarray): Weighted Shannon Entropy for all sequence indexes in the MSA
		"""

		## Recover the node using its ID from PhyloTree
		node:PhyloTree.__Node = self.tree.nodes[node_id]

		accessions:list = node.accessions

		n_j:int = len(accessions)
		log_n_j = np.log(n_j)

		## Occurrences of residues at position of i
		residue_counts:pd.DataFrame = self.msa.get_residue_counts(accessions)

		## Occurrences of AA at position i in the subfamily, minus gaps, normalized by number of proteins in subfamily
		## Array is of size number_of_residues x 20 (number of amino acids)
		a_i:pd.DataFrame = residue_counts.drop(columns='-')/n_j

		## ln(0) ==> bad things... instead make ln(1) ==> 0
		a_i[a_i==0] = 1

		se_ai:pd.series = -np.sum(a_i*np.log(a_i),axis=1)

		## Number of gaps at position in the subfamily
		G_i:pd.Series = residue_counts['-']

		se_Gi = G_i*(1/n_j)*np.log(1/n_j)

		## Normalization factor
		W_j:float = 1/log_n_j if n_j < 20 or not weigh else 1/np.log(20)

		## Normalize the combination of the AA and gap entropy
		E_i:pd.Series = W_j*(se_ai-se_Gi)

		return E_i

	def __calculate_average_entropy(self) -> np.array:

		nodal_entropy = self.__calculate_nodal_entropy()
		nodes_by_distance = self.tree.nodes_by_distance

		average_entropy:np.array = np.zeros(len(self.msa.msa.index))

		for idx,distance in enumerate(sorted(nodes_by_distance.keys())):

			subfamily_entropy:np.array = np.zeros(len(self.msa.msa.index))
			
			for node in nodes_by_distance[distance]:

				subfamily_entropy += nodal_entropy[node]

			average_entropy += 1/len(nodes_by_distance[distance])*subfamily_entropy

		average_entropy *= 1/len(nodes_by_distance.keys())

		return average_entropy

	def __calculate_nodal_entropy(self) -> dict:

		nodal_entropy:dict = {}
		tree_nodes:list = list(self.tree.nodes.keys())
		node_id:int

		for node_id in tree_nodes:

			if len(self.tree.nodes[node_id].accessions) == 1:

				nodal_entropy[node_id] = 0
				continue

			nodal_entropy[node_id] = self.__calculate_shannon_entropy(node_id)

		return nodal_entropy
	
	# def summarize_TEA_results(self,reference_accession:str|None=None) -> None:

	# 	if reference_accession is None:
	# 		reference_accession = self.msa.msa.columns[0]

	# 	outfile:pathlib.Path = self.outdir / "two_entropy.summary"

	# 	reference_seq_to_msa_indexes:np.array = self.msa.get_reference_indices(reference_accession)
	# 	consensus_sequence:pd = self.msa.get_residue_counts(self.msa.msa.columns)

	# 	with outfile.open('w') as OUT:

	# 		OUT.write(f"## Sequence_Pos\tMSA_Pos\tResidue\tGlobal_Shannon_Entropy\tAverage_Shannon_Entropy\tSequences_at_pos\tNum_of_gaps\n")

	# 		for res_index,msa_index in enumerate(reference_seq_to_msa_indexes):
				
	# 			global_entropy = self.global_entropy[msa_index]
	# 			average_entropy = self.average_entropy[msa_index]
	# 			num_of_gaps = consensus_sequence['-'][msa_index]
	# 			non_gapped_res = len(self.msa.msa.columns)

	# 			## Residue index
	# 			OUT.write(f"{res_index}")
	# 			## MSA index
	# 			OUT.write(f"\t{msa_index}")
	# 			## Global entropy
	# 			OUT.write(f"\t{global_entropy}")
	# 			## Average entropy
	# 			OUT.write(f"\t{average_entropy}")
	# 			## Residue of reference
	# 			# OUT.write(f"\t{SequenceUtilities.AAs[np.where(msa.loc[reference_id,'array'].T[msa_index]==1)[0][0]]}")
	# 			## Consenus
	# 			# residue_count = ";".join([f"{residue}:{int(count)}" for count,residue in sorted(list(zip(np.sum(msa.loc[msa.index,'array'],axis=0).T[msa_index],SequenceUtilities.AAs)),key=lambda x: x[0],reverse=True)])
	# 			# OUT.write(f"\t{residue_count}")
	# 			## Non-gapped residues
	# 			OUT.write(f"\t{non_gapped_res}")
	# 			## Gapped residues
	# 			OUT.write(f"\t{num_of_gaps}\n")

	# 	return None

	@property
	def global_entropy(self):
		if self.__global_entropy is None:
			self.__global_entropy = self.__calculate_shannon_entropy(0,weigh=False)
		return self.__global_entropy
	
	@property
	def average_entropy(self):
		if self.__average_entropy is None:
			self.__average_entropy = self.__calculate_average_entropy()
		return self.__average_entropy