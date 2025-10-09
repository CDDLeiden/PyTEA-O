#!/usr/bin/env python3

import numpy as np
import pandas as pd
from time import time
import typing
import heapq
import pathlib

from PyTEAO.utils.msa import MSA
from PyTEAO.trees.treebase import Tree
from PyTEAO.utils.general import valid_directory

class TwoEntropyAnalysis:

	def __init__(self,msa:MSA,tree:Tree,threads:int=1,outdir:str|pathlib.Path=pathlib.Path("./TEA").resolve(strict=False)):

		self.msa:MSA = msa
		self.residue_numbers = msa.residue_indexes
		self.tree:Tree = tree
		self.threads = threads
		self.outdir = valid_directory(outdir)
		self.__average_entropy = None
		self.__global_entropy = None
		self.__conservation_scores = None
		self.__specificity_scores = None

	def __calculate_shannon_entropy(self,node_id:int,weigh:bool=True):

		"""
		Calculates the weighted Shannon Entropy for the provided MSA utilizing the proposed equation by Ye, 2008

		Parameters
		msa (pd.DataFrame): DataFrame containing all residue counts for all sequence indexes for all sequences in the MSA

		Output
		E_i (np.ndarray): Weighted Shannon Entropy for all sequence indexes in the MSA
		"""

		## Recover the node using its ID from PhyloTree
		node:Tree.__Node = self.tree._nodes[node_id]

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

		average_entropy:np.array = np.zeros(len(self.residue_numbers))

		for idx,distance in enumerate(sorted(nodes_by_distance.keys())):

			subfamily_entropy:np.array = np.zeros(len(self.residue_numbers))
			
			for node in nodes_by_distance[distance]:

				subfamily_entropy += nodal_entropy[node]

			average_entropy += 1/len(nodes_by_distance[distance])*subfamily_entropy

		average_entropy *= 1/len(nodes_by_distance.keys())

		return average_entropy

	def __calculate_nodal_entropy(self) -> dict:

		nodal_entropy:dict = {}
		tree_nodes:list = list(self.tree._nodes.keys())
		node_id:int

		for node_id in tree_nodes:

			if len(self.tree._nodes[node_id].accessions) == 1:

				nodal_entropy[node_id] = 0
				continue

			nodal_entropy[node_id] = self.__calculate_shannon_entropy(node_id)

		return nodal_entropy

	def __calculate_conservation_scores(self):

		return -np.sqrt(np.square(self.average_entropy)+np.square(self.global_entropy))+0.5
	
	def __calculate_specificity_scores(self):

		return -np.sqrt(np.square(self.average_entropy)+np.square(self.global_entropy-1))+0.5

	@property
	def global_entropy(self):
		if self.__global_entropy is None:
			self.__global_entropy = self.__calculate_shannon_entropy(self.tree.root.node_id,weigh=False)
		return self.__global_entropy
	
	@property
	def average_entropy(self):
		if self.__average_entropy is None:
			self.__average_entropy = self.__calculate_average_entropy()
		return self.__average_entropy
	
	@property
	def conservation_scores(self):
		if self.__conservation_scores is None:
			self.__conservation_scores = self.__calculate_conservation_scores()
		return self.__conservation_scores

	@property
	def specificity_scores(self):
		if self.__specificity_scores is None:
			self.__specificity_scores = self.__calculate_specificity_scores()
		return self.__specificity_scores