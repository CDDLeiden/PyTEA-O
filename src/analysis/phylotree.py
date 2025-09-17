import numpy as np
import heapq
import typing

from src.utils.msa import MSA

class PhyloTree:

	class __Node:

		def __init__(self,node_id:int,accessions:list,branch_length:float):

			self.node_id:int = node_id
			self.accessions:list = accessions
			self.branch_length:float = branch_length
			self.children:list = []
			self.parent:PhyloTree.__Node = None


	def __init__(self,msa:MSA,threads:int=1) -> None:

		self.nodes:dict = {}
		self.distance_matrix = msa.distance_matrix
		self.next_node_id = 2*self.distance_matrix.shape[0]-2
		self.__root = None
		self.threads = threads

	def __add_node(self,accessions,branch_length) -> __Node:

		node = self.nodes[self.next_node_id] = self.__Node(
											node_id=self.next_node_id,
											accessions=accessions,
											branch_length=branch_length
										)
		
		self.next_node_id -= 1

		return node

	def __join_nodes(self,node_1,node_2,branch_length) -> __Node:

		new_node = self.__add_node(
			accessions=node_1.accessions+node_2.accessions,
			branch_length=branch_length
		)
		
		new_node.children = [node_1,node_2]

		node_1.parent = new_node
		node_2.parent = new_node

		return new_node
	
	def __pair_key(self,node_a_id:int,node_b_id:int) -> frozenset:

		return frozenset((node_a_id,node_b_id))
	
	def __init_upgma(self) -> typing.Tuple[dict,dict,list]:

		accessions:list = list(self.distance_matrix.columns)

		## Nodes that are available to be merged
		active_nodes:dict[int,PhyloTree.__Node] = {}
		
		## Number of accessions for a given node to weigh branch length contribution
		accession_to_nodeid:dict[str,int] = {}

		for accession in accessions:
			leaf:PhyloTree.__Node = self.__add_node(accessions=[accession],branch_length=0)
			active_nodes[leaf.node_id] = leaf
			accession_to_nodeid[accession] = leaf.node_id

		## Track the pairwise diatance between two accessions
		## Frozenset is used because it is immutable, therefore it can be set as a hash key, but also because sets
		## can be compared without needing to be ordered, whereas tuples need to be ordered
		## i.e. set(1,2) == set(2,1), tuple(1,2) != tuple(2,1)
		pairwise_distances:dict[frozenset,float] = {}

		## Using heap to avoid matrix manipulation/creation, O(n^2logn) vs O(n^3) complexity (faster)
		heap:list[tuple[float,int,int]] = []

		## Iterate over all accessions
		for i,a in enumerate(accessions):
			## Iterate over upper triangular
			for j in range(i+1,len(accessions)):
				b = accessions[j]
				ai = accession_to_nodeid[a]
				bi = accession_to_nodeid[b]
				d = float(self.distance_matrix.loc[a,b])
				pairwise_distances[self.__pair_key(ai,bi)] = d
				heapq.heappush(heap,(d,ai,bi))

		return active_nodes,pairwise_distances,heap

	def __upgma(self) -> __Node:

		active_nodes,pairwise_distances,heap = self.__init_upgma()
		
		while(len(active_nodes)>1):

			## Find the next active minimum distance
			while True:

				pairwise_distance,node_a_id,node_b_id = heapq.heappop(heap)
				
				## If both nodes are still active, proceed
				if node_a_id in active_nodes and node_b_id in active_nodes:
					break

			node_a:PhyloTree.__Node = active_nodes[node_a_id]
			node_b:PhyloTree.__Node = active_nodes[node_b_id]

			## Create new parent node
			new_node = self.__join_nodes(node_a,node_b,pairwise_distance/2)
			new_node_id = new_node.node_id

			## Remove no-longer active nodes
			del(active_nodes[node_a_id])
			del(active_nodes[node_b_id])

			## Activate new node
			active_nodes[new_node.node_id] = new_node

			node_a_weight = len(node_a.accessions)
			node_b_weight = len(node_b.accessions)

			for id in list(active_nodes.keys()):
				
				if id == new_node_id:
					continue
				
				node_a_dist = pairwise_distances.get(self.__pair_key(node_a_id,id))
				node_b_dist = pairwise_distances.get(self.__pair_key(node_b_id,id))

				new_node_dist = (node_a_dist*node_a_weight+node_b_dist*node_b_weight)/(node_a_weight+node_b_weight)
				
				## Add new pairwise distance
				pairwise_distances[self.__pair_key(new_node_id,id)] = new_node_dist

				## Add new distance to the heap
				heapq.heappush(heap,(new_node_dist,new_node_id,id))

		return next(iter(active_nodes.values()))

	@property
	def root(self):
		if self.__root is None:
			self.__root = self.__upgma()
		return self.__root