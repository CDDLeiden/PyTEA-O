import numpy as np
import heapq
import typing

from src.utils.msa import MSA

class PhyloTree:

	class __Node:

		def __init__(self,node_id:int,accessions:list,node_placement:float):

			self.node_id:int = node_id
			self.accessions:list = accessions
			self.node_placement:float = node_placement
			self.children:list[int] = []
			self.parent:PhyloTree.__Node = None


	def __init__(self,msa:MSA,threads:int=1) -> None:

		self.nodes:dict = {}
		self.distance_matrix = msa.distance_matrix
		self.next_node_id = 2*self.distance_matrix.shape[0]-2
		self.threads = threads
		self.root = self.__upgma()
		
		self.__nodes_by_distance = None

	def __add_node(self,accessions,node_placement) -> __Node:

		node = self.nodes[self.next_node_id] = self.__Node(
											node_id=self.next_node_id,
											accessions=accessions,
											node_placement=node_placement
										)
		
		self.next_node_id -= 1

		return node

	def __join_nodes(self,node_1:__Node,node_2:__Node,node_placement) -> __Node:

		new_node = self.__add_node(
			accessions=node_1.accessions+node_2.accessions,
			node_placement=node_placement
		)
		
		new_node.children = [node_1.node_id,node_2.node_id]

		node_1.parent = new_node.node_id
		node_2.parent = new_node.node_id

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
			leaf:PhyloTree.__Node = self.__add_node(accessions=[accession],node_placement=0)
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
	
	def __get_nodes_by_distance(self) -> dict:

		root:PhyloTree.__Node = self.root
		root_id:int = root.node_id
		root_placement:float = root.node_placement
		
		num_of_leafs:int = len(self.distance_matrix.columns)
		
		distances:dict = {}
		active_nodes:dict[int:PhyloTree.__Node] = {root_id:root}

		heap = []
		heapq.heappush(heap,(0,root_id))

		parent_node_placement:float
		parent_node_id:int

		while len(active_nodes) < num_of_leafs-1:

			parent_node_placement,parent_node_id = heapq.heappop(heap)

			parent_node:PhyloTree.__Node = self.nodes[parent_node_id]

			child_node_id:int

			for child_node_id in parent_node.children:

				child_node:PhyloTree.__Node = self.nodes[child_node_id]

				## Traversing tree in reverse, need to invert node_placement by root_placement
				child_node_placement = root_placement - child_node.node_placement

				## Add the new child nodes to active nodes
				active_nodes[child_node_id] = child_node

				heapq.heappush(heap,(child_node_placement,child_node_id))

			## Parent node no longer present
			del(active_nodes[parent_node_id])

			## Track which nodes are at current placement
			distances[abs(parent_node_placement-root_placement)] = list(active_nodes.keys())

		return distances
	
	@property
	def nodes_by_distance(self):
		if self.__nodes_by_distance is None:
			self.__nodes_by_distance = self.__get_nodes_by_distance()
		return self.__nodes_by_distance