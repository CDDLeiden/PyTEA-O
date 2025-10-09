import heapq
import typing
import pathlib

from PyTEAO.trees.treebase import Tree
from PyTEAO.utils.msa import MSA

class PhyloTree(Tree):

	def __init__(self,msa:MSA) -> None:

		super().__init__(msa)

		self.distance_matrix = msa.distance_matrix

	def define_tree(self,lineage_file) -> None:

		lineage_file = pathlib.Path(lineage_file)

		self._nodes = {}

		max_placement = 0

		with lineage_file.open('r') as IN:

			line_num = 0
			
			for line in IN:

				line_num += 1

				line = line.strip()

				if line[0] == "#":
					continue
				
				node_id,node_placement,accessions = line.split("\t")

				accessions = accessions.split(",")
				node_id,node_placement = [int(node_id),int(node_placement)]

				node = self.define_node(accessions,node_placement,node_id)

				if self._nodes.get(node_id) is not None:
					raise ValueError(f"Duplicate node_id [{node_id}] passed on line {line_num}\n\t{line}\nTree building halted to prevent SNAFUs.")
				
				self._nodes[node.node_id] = node

				max_placement = max_placement if node_placement < max_placement else node_placement

		self._root = self.join_nodes(list(self._nodes.values()),max_placement+1)
	
	def __pair_key(self,node_a_id:int,node_b_id:int) -> frozenset:

		return frozenset((node_a_id,node_b_id))
	
	def __init_upgma(self) -> typing.Tuple[dict,dict,list]:

		accessions:list = list(self.distance_matrix.columns)

		active_nodes,accession_to_node_id = self.bulk_node_creation(accessions)

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
				ai = accession_to_node_id[a]
				bi = accession_to_node_id[b]
				d = float(self.distance_matrix.loc[a,b])
				pairwise_distances[self.__pair_key(ai,bi)] = d
				heapq.heappush(heap,(d,ai,bi))

		return active_nodes,pairwise_distances,heap

	def __upgma(self) -> Tree.Node:

		active_nodes,pairwise_distances,heap = self.__init_upgma()
		
		while(len(active_nodes)>1):

			## Find the next active minimum distance
			while True:

				pairwise_distance,node_a_id,node_b_id = heapq.heappop(heap)
				
				## If both nodes are still active, proceed
				if node_a_id in active_nodes and node_b_id in active_nodes:
					break

			node_a:Tree.Node = active_nodes[node_a_id]
			node_b:Tree.Node = active_nodes[node_b_id]

			## Create new parent node
			new_node = self.join_nodes([node_a,node_b],pairwise_distance/2)
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
		if self._root is None:
			self._root = self.__upgma()
		return self._root