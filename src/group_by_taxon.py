#!/usr/bin/env python3

from os.path import isfile, abspath
import os
from subprocess import run
from time import time
from asyncio.queues import Queue
from asyncio import create_task, gather, run
from typing import Generator
import mmap


def setup_db(db_dir:str=None) -> dict:

	"""
	Setup the necessary databases for protein accession number to taxonoic IDs,
	downloading files as necessary.

	Args
		db_dir	Directory containing necessary datases

	Return
		database_files	Dictionary object linking to database file locations
	"""

	pipeline_path = "/".join(abspath(__file__).split("/")[0:-2])
	db_dir = f"{pipeline_path}/DB_files"

	database_files = {
		'names':f"{db_dir}/names.dmp",
		'nodes':f"{db_dir}/nodes.dmp",
		'a2tax':f"{db_dir}/prot.accession2taxid.FULL.gz"
	}

	if not isfile(database_files['names']) and not isfile(database_files['nodes']):

		###########################################################################################
		# Download TaxDump File
		###########################################################################################

		run(['wget','https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz','-O',f"{db_dir}/taxdump.tar.gz"])
		run(['tar','-xzvf',f"{db_dir}/taxdump.tar.gz","--directory",f"{db_dir}"])
		excess_dmps = [
			'taxdump.tar.gz',
			'citations.dmp',
			'delnodes.dmp',
			'division.dmp',
			'gc.prt',
			'gencode.dmp',
			'images.dmp',
			'merged.dmp',
			'readme.txt'
		]
		for item in excess_dmps:
			run(['rm',f"{db_dir}/{item}"])


	if not isfile(database_files['a2tax']):
		
		###########################################################################################
		# Download Accession2Taxid File (~18Gb)
		###########################################################################################

		ftp_link = f"https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.FULL.gz"
		connections = 10
		run([
			'aria2c',
			'-s',f"{connections}",
			'-x',f"{connections}",
			'-d',f"{db_dir}",
			'--continue',
			f"{ftp_link}"
		])

		###########################################################################################
		# Download Accession2Taxid MD5 File
		###########################################################################################

		run([
			'aria2c',
			'-s',f"{connections}",
			'-x',f"{connections}",
			'-d',f"{db_dir}",
			'--continue',
			f"{ftp_link}.md5"
		])

	# checksum(f"{db_dir}/{accession_file}.md5")

	return database_files

def load_nodes(nodes_file:str=None) -> dict:

	"""
	Parses taxonomic database file, linking child IDs to parent IDs

	Args
		nodes_file	Path to NCBI nodes.dmp database file

	Returns
		node_links	Dictionary object linking child taxonomic ID to parent taxonmic ID
	"""

	node_links = {}

	with open(nodes_file,'r') as IN:

		for line in IN:

			line = line.strip()

			line_data = line.split("\t|\t")

			child_org_tid,parent_org_tid,child_rank = line_data[0:3]

			node_links[int(child_org_tid)] = {'rank':child_rank,'parent_taxid':int(parent_org_tid)}

	return node_links

def get_sequence_accessions(msa_file:str=None) -> set:

	"""
	Extracts accession numbers from Multiple Sequence Alignment (MSA) files
	in FASTA format.

	Args
		msa_file	Multiple Sequence Alignment (MSA) file in FASTA format

	Returns
		accessions_set	Set of all accession numbers in MSA file
	"""

	accessions_set = set()

	with open(msa_file,'r') as IN:
		for line in IN:
			line = line.strip()
			if line == "":
				next
			if line[0] == ">":
				accessions_set.add(line[1:])

	return accessions_set

async def read_a2tax_file(a2tax_file:str=None,set_chunk_size:int=(5*(1024**2)),input_queue:Queue=None)-> None:

	from gzip import open as gopen

	## Get size of accessions_to_taxid file
	file_size = os.stat(a2tax_file).st_size

	chunk_index = 0

	## If file is smaller than chunk size, set chunks size to file size
	set_chunk_size = set_chunk_size if file_size > set_chunk_size else file_size

	## Create binary file handle to read from
	with gopen(a2tax_file,'rt') as IN:

		## Create memorymap object to copy from, accessing the full file
		with mmap.mmap(IN.fileno(),length=0,access=mmap.ACCESS_COPY) as MM_IN:
			
			## Iterate over each chunk in the file
			for file_pointer_start,read_chunk_size in read_file_chunks(file_handle=MM_IN,chunk_size=set_chunk_size):
				
				chunk_index += 1

				if chunk_index%10 == 0:
					
					print(f"Processed {chunk_index} items ({chunk_index*set_chunk_size/(1024**3):.2f})")
				
				## Move file pointer to chunk start
				MM_IN.seek(file_pointer_start)

				## Grab file chunk
				file_chunk = IN.read(read_chunk_size)

				start = time()

				## Signal EOF and breakout of the generator
				if file_chunk is None:
					
					await input_queue.put(None)
					break

				if file_chunk == "":
					
					await input_queue.put(None)
					break

				if file_chunk[-1] != "\n":
					file_chunk += IN.readline()

				## Wait for the input queue to have an available spot
				await input_queue.put(file_chunk)

	await input_queue.put(None)
	return None

def read_file_chunks(file_handle:mmap.mmap=None,chunk_size:int=(5*(1024**2))) -> Generator[None,int,int]:

	"""

	Args
		file_handle	Memory map (mmap.mmap) file object
		chunk_size	Size of chunks to create in byets [Default: 5*(1024**2) bytes]

	Returns
		file_point_start	Address of current position in memory map object
		file_chunk_length	Length of memory map object chunk

	"""

	while True:

		## Grab file pointer position
		file_pointer_start = file_handle.tell()

		## Move file pointer position to next chunk
		file_handle.seek(chunk_size,1)

		## Return file start address and file chunk size
		yield file_pointer_start, file_handle.tell() - file_pointer_start

async def convert(file_chunk=None,accession_set:set=None,a2tax_links:dict=None):

	if file_chunk is None:
		return None
	
	if len(accession_set) == 0:
		return a2tax_links

	for line in file_chunk.split("\n"):
		if line == "":
			continue
		accession,taxid = line.strip().split("\t")
		# print(f"Accession: {accession} | TaxID: {taxid}"
		if accession in accession_set:
			try:
				a2tax_links[accession] = taxid
				accession_set.discard(accession)
			except:
				pass

	return a2tax_links

async def convert_accession_to_taxid(
			input_queue:Queue=None,
			a2tax_links:dict=None,
			accession_set:set=None,
			outdir:str="TaxonomyGrouping",
	) -> None:

	while True:
		
		file_chunk = await input_queue.get()

		input_queue.task_done()

		if file_chunk is not None:

			await convert(
					file_chunk=file_chunk,
					accession_set=accession_set,
					a2tax_links=a2tax_links
			)

		else:
			
			break

	with open(f"{outdir}/accession_taxids.tsv",'w') as OUT:

		for key,value in a2tax_links.items():
			
			OUT.write(f"{key}\t{value}\n")

	return None

async def link_accessions_to_taxids(
			a2tax_file:str=None,
			accession_set:set=None,
			max_chunks:int=10,
			chunk_size:int=(5*(1024**2))
	) -> dict:

	"""
	Generates links between protein accessions and taxonomy IDs using NCBIs prot.accession2taxid.FULL.gz file.

	|Args
	|	a2tax_file: Path to NCBI's prot.accession2taxid.FULL.gz file
	|	accession_set: Unique set of protein accessions acquired from MSA file
	|	max_chunk: Number of file chunks to load in at a single time
	|	chunk_size: Size of file chunks
	"""

	a2tax_links = {}

	## Create chunk queue
	chunk_queue = Queue(maxsize=max_chunks)

	## Create chunk reading task
	read_a2tax = create_task(
					read_a2tax_file(
						a2tax_file=a2tax_file,
						set_chunk_size=chunk_size,
						input_queue=chunk_queue
					)
	)

	## Create chunk processing task
	conversion_task = create_task(
						convert_accession_to_taxid(
							input_queue=chunk_queue,
							a2tax_links=a2tax_links,
							accession_set=accession_set
						)
	)

	await gather(read_a2tax,conversion_task)

	return a2tax_links

def group_accessions(
		a2tax_links:dict=None,
		nodes_links:dict=None,
		outdir:str="TaxonomyGrouping",
	) -> None:

	from time import time

	database_files = setup_db()
	nodes_links = load_nodes(nodes_file=database_files['nodes'])
	accession_set = get_sequence_accessions(msa_file=args.msa)

	start = time()

	a2tax_links = {}

	run(
		link_accessions_to_taxids(
			a2tax_file=database_files['a2tax'],
			accession_set=accession_set,
			max_chunks=args.num_chunks,
			chunk_size=args.chunk_size*(1024**2)
		)
	)

	# with open("test.txt",'r') as IN:
	# 	for line in IN:
	# 		accession,taxid = line.strip().split("\t")
	# 		a2tax_links[accession] = int(taxid)

	taxonomy = {}

	for accession in a2tax_links.keys():

		taxid = a2tax_links[accession]
		
		while taxid:

			try:
				parent_tid = nodes_links[taxid]['parent_taxid']
				rank = nodes_links[taxid]['rank']
			except:
				break

			if rank not in taxonomy.keys():
				taxonomy[rank] = {}

			if taxid not in taxonomy[rank].keys():
				taxonomy[rank][taxid] = []

			taxonomy[rank][taxid].append(accession)

			taxid = parent_tid

			if taxid == 1:
				break

	desired_outputs = [
		"superkingdom",
		"kingdom",
		"phylum",
		"class",
		"order",
		"family",
		"genus",
		"species"
	]

	groupings = {}
	counts = {}

	for rank in taxonomy.keys():

		counts[rank] = 0

		for index,taxid in enumerate(taxonomy[rank].keys()):

			for accession in taxonomy[rank][taxid]:

				if accession not in groupings.keys():

					groupings[accession] = {}

				groupings[accession][rank] = hex(index)

				counts[rank] += 1

	with open(f"{outdir}/groupings.csv",'w') as OUT:
		
		headers = ";".join(sorted(groupings.keys()))

		OUT.write(f"## {headers}\n")

		for rank in desired_outputs:

			line = ";".join([groupings[x][rank] if rank in groupings[x].keys() else "N/A" for x in sorted(groupings.keys())])

			OUT.write(f"{line}\n")

	return None

def main(args:dict=None) -> int:

	from time import time

	start = time()

	## Setup database files
	database_files = setup_db()

	## Link child taxids to parent taxids
	nodes_links = load_nodes(
		nodes_file=database_files['nodes']
	)
	
	## Extract MSA accession numbers
	accession_set = get_sequence_accessions(
		msa_file=args.msa
	)

	## Link MSA accession numbers to taxids
	a2tax_links = run(
		link_accessions_to_taxids(
			a2tax_file=database_files['a2tax'],
			accession_set=accession_set,
			max_chunks=args.num_chunks,
			chunk_size=args.chunk_size*(1024**2)
		)
	)

	## Group accession numbers by taxids
	group_accessions(
		a2tax_links=a2tax_links,
		nodes_links=nodes_links
	)

	return None

if __name__ == "__main__":

	from argparse import ArgumentParser

	GetOptions = ArgumentParser()

	GetOptions.add_argument(
					"-m","--msa",
					required=True,
					help="Path to MSA-file containing accessions to group by taxonomy."
	)

	GetOptions.add_argument(
					"-s","--chunk_size",
					required=False,
					type=int,
					default=5,
					help="Size of file chunks in megabases (Mb) [Default: 5Mb]"
	)
	
	GetOptions.add_argument(
					"-n","--num_chunks",
					required=False,
					type=int,
					default=10
	)

	GetOptions.add_argument(
					"-t","--threads",
					required=False,
					type=int,
					default=2
	)

	GetOptions.add_argument(
					"-o","--outdir",
					type=str,
					default="TaxonomyGrouping"
	)

	main(args=GetOptions.parse_known_args()[0])