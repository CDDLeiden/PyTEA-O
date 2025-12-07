import pathlib

from PyTEAO.trees.treebase import Tree

class UserTree(Tree):

	def __init__(self,msa,lineage_file):

		super().__init__(msa)

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

		self.root = self.join_nodes(list(self._nodes.values()),max_placement+1)