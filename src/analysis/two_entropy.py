import pandas as pd
import numpy as np
import time
import os
import argparse
import pathlib

from src.utils.general import read_subfamilies, create_directory, merge_commflags_with_kwargs, valid_file, parse_args, valid_directory,get_file_name
from src.utils.msa import load_msa, get_reference_indices, get_consensus_sequence
from src.analysis.shannon_entropy import shannon_entropy_multiprocess,superfamily_shannon_entropy, calculate_subfamily_shannon_entropy
from src.grouping.upgma_subfam_grouping import generate_subgroups
from src.utils.sequence import SequenceUtilities
from src.analysis.upgma import run as run_upgma


ARGUEMENTS = {
	'm':{
		'flag':'msa_file',
		'required':True,
		'help':"File path to Multiple Sequence Alignment",
		'type':valid_file,
		'CLI': True,
		'call': True
	},
	's': {
		'flag': 'subfamilies_file',
		'help': "File location of subfamily assignment (TXT file)",
		'type': valid_file,
		'CLI': True,
		'call': True
	},
	'o': {
		'flag': 'outdir',
		'help': "Output file location",
		'type': str,
		'default': "TEA",
		'type': valid_directory,
		'CLI': True,
		'call': True
	},
	'p': {
		'flag': 'program_mode',
		'help': "Run program in standard TEA or objective TEAO",
		'type': str,
		'choices': ['TEA', 'TEAO'],
		'default': "TEAO",
		'CLI': True,
		'call': True
	},
	'r': {
		'flag': 'reference_id',
		'help': "Reference ID(s)",
		'type': str,
		'CLI': True,
		'call': True
	},
	'x': {
		'flag': 'similarity_matrix',
		'help': "Calculate scores according to provided matrix",
		'type': str,
		'CLI': True,
		'call': True
	},
	't': {
		'flag': 'threads',
		'help': "Number of threads",
		'type': int,
		'default': 1,
		'CLI': True,
		'call': True
	}
}

def read_temp_tea_file(temp_file:str) -> dict:

	results = {}

	with open(temp_file,'r') as IN:
		for line in IN:
			line = line.strip()
			if line == "":
				continue
			key,data = line.split("\t")
			results[int(key)] = [float(x) for x in data.split(",")]

	return results

def summarize_tea_results(reference_id:str,msa:pd.DataFrame,global_SE:pd.DataFrame,average_SE:np.ndarray,residues:np.ndarray,gaps:np.ndarray,outdir:str) -> None:

	print(f"\n\tSummarizing Two Entropy Analysis for reference sequence...")

	if not reference_id: 
		print(f"\n\t\t[N] Reference ID not provided, defaulting to first sequence in MSA [{msa.index[0]}]...")
		reference_id = msa.index[0]

	outfile = outdir / "two_entropy.summary"

	reference_seq_to_msa_indexes = get_reference_indices(ref=reference_id,msa=msa)

	get_consensus_sequence(msa)

	with outfile.open('w') as OUT:
		OUT.write(f"## Sequence_Pos\tMSA_Pos\tResidue\tGlobal_Shannon_Entropy\tAverage_Shannon_Entropy\tSequences_at_pos\tNum_of_gaps\n")

		for res_index,msa_index in enumerate(reference_seq_to_msa_indexes):
			
			global_entropy = global_SE[0][msa_index]
			average_entropy = average_SE[0][msa_index]
			non_gapped_res = int(np.sum(residues[:,msa_index]))
			num_of_gaps = int(gaps[0][msa_index])

			## Residue index
			OUT.write(f"{res_index}")
			## MSA index
			OUT.write(f"\t{msa_index}")
			## Global entropy
			OUT.write(f"\t{global_entropy}")
			## Average entropy
			OUT.write(f"\t{average_entropy}")
			## Residue of reference
			OUT.write(f"\t{SequenceUtilities.AAs[np.where(msa.loc[reference_id,'array'].T[msa_index]==1)[0][0]]}")
			## Consenus
			residue_count = ";".join([f"{residue}:{int(count)}" for count,residue in sorted(list(zip(np.sum(msa.loc[msa.index,'array'],axis=0).T[msa_index],SequenceUtilities.AAs)),key=lambda x: x[0],reverse=True)])
			OUT.write(f"\t{residue_count}")
			## Non-gapped residues
			OUT.write(f"\t{non_gapped_res}")
			## Gapped residues
			OUT.write(f"\t{num_of_gaps}\n")

	return None

def TEA(subfamilies_file:str,outdir:pathlib.Path,msa:pd.DataFrame,threads:int,mode:str,msa_file:str,file_prefix:str="TEA") -> np.ndarray:


	subfamilies:dict = read_subfamilies(subfamilies_file)

	(outdir := outdir/'.data').mkdir(mode=0o755,parents=True,exist_ok=True)

	num_of_MSA_res,_ = msa.shape

	for file in [x for x in outdir.iterdir() if x.is_file()]:
		subfamilies.pop(get_file_name(file).replace("_","."),None)

	# if temp_entropy_file.is_file():
	# 	results = read_temp_tea_file(temp_entropy_file)

	# shannon_entropy_multiprocess(msa=msa,temp_entropy_file=temp_entropy_file,results=results,subfamilies=subfamilies,keys=keys,threads=threads)

	# results = read_temp_tea_file(temp_entropy_file)

	# E_i = np.zeros(num_of_MSA_res)
	
	# for index in results.keys():
	# 	E_i += results[index]
	# E_i /= len(results.keys())

	# with open(f"{outdir}/{file_prefix}.average_entropy",'w') as OUT:
	# 	for msa_pos,se in enumerate(E_i):
	# 		OUT.write(f"{msa_pos}\t{se:.4f}\n")

	# return E_i

	return

def run(args:argparse.Namespace|None=None,**kwargs) -> pathlib.Path:

	"""
	Performs the (T)wo (E)ntropy (A)nalysis analysis on a given MSA, finding the global entropy,
	and the average entropy across subfamilies.

	Parameters:
	msa_file (str): File containing MSA for protein of interest
	subfamily_file (str): File containing user defined subfamily groupings for proteins in the provided MSA
	
	Returns:
	pathlib.Path pointing to summarized TEA results
	"""

	print(f"\n\tTwo Entropy Analysis started")



	args:argparse.Namespace = merge_commflags_with_kwargs(
		cli_args=args,
		function_args=ARGUEMENTS,
		**kwargs
	)

	(outdir := args.outdir/'TEA').mkdir(mode=0o755,parents=True,exist_ok=True)

	msa:pd.DataFrame = load_msa(
		args.msa_file,
		outdir=outdir
	)
	
	gE_i,a_i,G_i = superfamily_shannon_entropy(
		msa=msa,
		outdir=outdir
	)

	if (args.program_mode == "TEA" and not args.subfamilies_file) or args.program_mode == "TEAO":

		if args.program_mode == "TEA":
			print("\n\n\t[W] --program_mode specified TEA algorithm, but --subfamilies_file not provided. Defaulting to TEAO mode...\n")
		phylo_tree = run_upgma(msa_dataframe=msa,threads=args.threads,outdir=outdir)
		args.subfamilies_file = outdir/'UPGMA'/'upgma.subfamilies'
		generate_subgroups(tree=phylo_tree,outfile=args.subfamilies_file)

	aE_i:np.array = TEA(
		subfamilies_file=args.subfamilies_file,
		outdir=outdir,
		msa=msa,
		threads=args.threads,
		mode=args.program_mode,
		msa_file=args.msa_file
	)

	results_file:pathlib.Path = summarize_tea_results(
		reference_id=args.reference_id,
		msa=msa,
		global_SE=gE_i,
		average_SE=aE_i,
		residues=a_i,
		gaps=G_i,
		outdir=args.outdir
	)

	return results_file

if __name__ == '__main__':

	run(args=parse_args(function_args=ARGUEMENTS))