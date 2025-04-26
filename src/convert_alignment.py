#!/usr/bin/env python3

from argparse import ArgumentParser
from os import makedirs
from textwrap import wrap

def read_MSA(msa_file:str=None) -> tuple[dict,str]:

	"""Reads an MSA from a given file, and converts it to a dictionary object, with loci and sequences as keys and
	values, repectively.

	Inputs:
		MSA_File: File path to MSA, either in FASTA or Clustal format.

	Returns:
		MSA: Dictionary object containig MSAs
	"""

	file_prefix:str = msa_file.split("/")[-1].split(".")[0]

	msa = {}
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

	msa = {x:"".join(msa[x]) for x in list(msa.keys())}

	return file_prefix,msa

def to_one_liner(aligned_sequences:dict,output:str,reference:str) -> None:

	with open(output,'w') as OUT:

		if reference:

			OUT.write(f">{reference}\n{aligned_sequences[reference]}\n")

		for sequence in sorted(aligned_sequences.keys()):

			if sequence == "reference":
				next

			OUT.write(f">{sequence}\n{aligned_sequences[sequence]}\n")

	return None

def extract(aligned_sequences:dict,outdir=str) -> None:

	fasta_dir = f"{outdir}/FASTAs" 
	makedirs(f"{outdir}/FASTAs",0o755,exist_ok=True)

	for locus in aligned_sequences.keys():

		with open(f"{fasta_dir}/{locus}.fasta",'w') as FASTA:

			sequence = "\n".join(wrap(aligned_sequences[locus].replace("-",""),60))

			FASTA.write(f">{locus}\n{sequence}\n")

	return None

def convert(aligned_sequences:dict,align_type:str,outdir:str,file_prefix:str,reference:str) -> None:

	if reference and reference not in aligned_sequences.keys():

		print(f"\n[W]  Reference locus - {reference} - was not found in the alignment set. Defaulting to normal functionality.\n")

		reference = None

	if align_type == "one_liner":
		
		to_one_liner(aligned_sequences=aligned_sequences,output=f"{outdir}/{file_prefix}.ola",reference=reference)
	
	if align_type == "clustal":
		...
	
	if align_type == "fasta":
		...

	return

def setup(outdir:str) -> None:
	
	makedirs(outdir,0o750,exist_ok=True)

	return None

def run(msa_file:str,outdir:str,extract_fastas:str,desired_alignment_type:str,reference:str) -> int:

	setup(outdir=outdir)

	file_prefix:str
	alignment_data:dict

	file_prefix,alignment_data = read_MSA(msa_file=msa_file)

	convert( aligned_sequences=alignment_data,
			 align_type=desired_alignment_type,
			 outdir=outdir,
			 file_prefix=file_prefix,
			 reference=reference
	)

	if extract_fastas: extract(aligned_sequences=alignment_data,outdir=outdir)

	return 0

if __name__ == "__main__":

	GetOptions = ArgumentParser()

	GetOptions.add_argument("-m","--msa_file",required=True,type=str)
	GetOptions.add_argument("-o","--outdir",default="alignment_conversion")
	GetOptions.add_argument("-x","--extract_fastas",default=False,action='store_true')
	GetOptions.add_argument("-a","--desired_alignment_type",default=["one_liner"],choices=["one_liner,clustal,fasta"])
	GetOptions.add_argument("-r","--reference",default=None,type=str)

	run(GetOptions.parse_known_args()[0])