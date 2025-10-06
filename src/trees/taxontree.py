import requests
from xml.etree import ElementTree
from tenacity import retry, wait_fixed, stop_after_attempt
import time
import warnings

from src.trees.treebase import Tree
from src.utils.msa import MSA

class TaxonTree(Tree):

	__BATCH_SIZE:int = 200
	__WAIT_TIME:float = 1/3
	__RANKS:list = ["domain","kingdom","phylum","class","order","family","genus","species"]

	def __init__(self,msa:MSA) -> None:

		super().__init__(msa)

		self.root = self.build_tree()

	@classmethod
	def protein_accessions_to_taxids(cls,accessions) -> dict:

		""" Converts provided list of protein sequence accessions into their respective NCBI taxids using
		NCBI Entrez Eutils. Accessions that do not link to NCBI taxids will return _None_.

		Returns:
			dict: Protein Accession -> NCBI TaxId
		"""

		@retry(wait=wait_fixed(cls.__WAIT_TIME), stop=stop_after_attempt(3))
		def __get_taxid(accessions):
			
			esummary_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
			params = {
				"db": "protein",
				"id": ",".join(accessions),
				"retmode": "xml"
			}
			
			r = requests.get(esummary_url, params=params)
			r.raise_for_status()

			return ElementTree.fromstring(r.content)

		def __parse_taxids(content):

			links:dict = {}

			for docsum in content.findall(".//DocSum"):
				accession = docsum.findtext(".//Item[@Name='Caption']")
				taxid = docsum.findtext(".//Item[@Name='TaxId']")
				links[accession] = taxid

			return links
		
		links:dict = {}

		for batch in [accessions[i:i+cls.__BATCH_SIZE] for i in range(0,len(accessions),cls.__BATCH_SIZE)]:

			content = __get_taxid(batch)

			links.update(__parse_taxids(content))

			time.sleep(cls.__WAIT_TIME)

		return links
	
	@classmethod
	def taxid_to_taxonomic_lineage(cls,taxids) -> dict:

		""" Expands provided list of valid NCBI TaxIds to their respective full taxonomic lineage. Invalid TaxIds will
		return _None_.

		Returns:
			dict: Taxid -> Full taxonomic lineage
		"""

		@retry(wait=wait_fixed(cls.__WAIT_TIME), stop=stop_after_attempt(3))
		def __get_lineage(taxids) -> ElementTree:
			
			esummary_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
			
			params = {
				"db": "taxonomy",
				"id": ",".join(taxids),
				"retmode": "xml"
			}
			
			r = requests.get(esummary_url, params=params)
			
			r.raise_for_status()

			return ElementTree.fromstring(r.content)

		def __parse_lineage(content:ElementTree):

			lineages:dict = {}

			for entry in content.findall("Taxon"):

				taxid = entry.findtext("TaxId")

				lineages[taxid] = {rank:None for rank in cls.__RANKS}

				lineage = entry.find("LineageEx")

				for taxon in lineage.findall("Taxon"):

					rank = taxon.findtext("Rank")

					if rank not in cls.__RANKS:
						continue

					lineages[taxid][rank] = taxon.findtext("TaxId")

			return lineages
		
		lineages = {}
		
		for batch in [taxids[i:i+cls.__BATCH_SIZE] for i in range(0,len(taxids),cls.__BATCH_SIZE)]:

			content = __get_lineage(batch)

			lineages.update(__parse_lineage(content))

			time.sleep(cls.__WAIT_TIME)

		for taxid in lineages.keys():

			lineages[taxid]['species'] = taxid
		
		return lineages
	
	@classmethod
	def protein_accessions_to_taxonomic_lineages(cls,accessions):
		
		links:dict = cls.protein_accessions_to_taxids(accessions)

		lineages:dict = cls.taxid_to_taxonomic_lineage(list(links.values()))

		return {x:lineages[links[x]] for x in links.keys()}

	def build_tree(self):

		accessions = self.msa.accessions

		lineages = self.protein_accessions_to_taxonomic_lineages(accessions)

		active_nodes,accessions_to_nodes_ids = self.bulk_node_creation(accessions)

		for ridx,rank in enumerate(reversed(self.__RANKS),start=1):

			taxids = {}

			node:Tree.Node

			for node in active_nodes.values():

				taxid = lineages[node.accessions[0]][rank]

				## Skip ranks that are not identified, don't falsely group into 'None' taxid
				if taxid is None:
					continue

				if taxids.get(taxid) is None:

					taxids[taxid] = []

				taxids[taxid].append(node)

			for nodes in taxids.values():

				## Don't create a new node for a solo taxid, maintain node in roster
				if len(nodes) < 2:
					continue

				## Merge previous nodes and add it to the active nodes
				new_node:TaxonTree.Node = self.join_nodes(nodes,ridx)
				active_nodes[new_node.node_id] = new_node

				## Remove merged nodes from active nodes
				for node in nodes:
					del(active_nodes[node.node_id])

		return next(iter(active_nodes.values()))