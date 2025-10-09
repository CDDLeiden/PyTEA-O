import heapq

class Tree:

	class Node:

		def __init__(self,node_id:int,accessions:list,node_placement:float):

			self.node_id:int = node_id
			self.accessions:list = accessions
			self.node_placement:float = node_placement
			self.children:list[int] = []
			self.parent:Tree.Node = None

	def __init__(self,msa):
		
		self.msa = msa
		self._nodes:dict = {}
		self.next_node_id:int = 0
		self._root = None
		self.__nodes_by_distance = None

	def create_node(self,accessions,node_placement,node_id) -> Node:

		node = self.Node(
						node_id=node_id,
						accessions=accessions,
						node_placement=node_placement
					)

		return node

	def add_node(self,accessions,node_placement) -> Node:

		node = self._nodes[self.next_node_id] = self.create_node(accessions,node_placement,self.next_node_id)

		self.next_node_id += 1

		return node
	
	def define_node(self,accessions,node_placement,node_id) -> Node:

		self.next_node_id = node_id + 1

		return self.create_node(accessions,node_placement,node_id)
	
	def bulk_node_creation(self,accessions) -> dict[str,int]:

		## Nodes that are available to be merged
		node_id_to_nodes:dict[int,Tree.Node] = {}
		
		## Number of accessions for a given node to weigh branch length contribution
		accession_to_node_id:dict[str,int] = {}

		for accession in accessions:
			leaf:Tree.Node = self.add_node(accessions=[accession],node_placement=0)
			node_id_to_nodes[leaf.node_id] = leaf
			accession_to_node_id[accession] = leaf.node_id

		return node_id_to_nodes,accession_to_node_id
	
	def join_nodes(self,nodes:list[Node],node_placement:int) -> Node:

		accessions = []

		for node in nodes:

			accessions += node.accessions
			node.parent = self.next_node_id

		new_node = self.add_node(
			accessions=accessions,
			node_placement=node_placement
		)

		new_node.children = [node.node_id for node in nodes]

		return new_node
	
	def get_nodes_by_distance(self) -> dict:

		root:Tree.Node = self.root
		root_id:int = root.node_id
		root_placement:float = root.node_placement
		
		num_of_leafs:int = len(root.accessions)
		
		distances:dict = {}
		active_nodes:dict[int:Tree.Node] = {root_id:root}

		heap = []
		heapq.heappush(heap,(0,root_id))

		parent_node_placement:float
		parent_node_id:int

		while len(heap) > 0:

			parent_node_placement,parent_node_id = heapq.heappop(heap)

			parent_node:Tree.Node = self._nodes[parent_node_id]

			child_node_id:int

			if len(parent_node.children) == 0:

				continue

			for child_node_id in parent_node.children:

				child_node:Tree.Node = self._nodes[child_node_id]

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
			self.__nodes_by_distance = self.get_nodes_by_distance()
		return self.__nodes_by_distance