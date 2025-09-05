import pandas as pd
import numpy as np
import time
import os
import argparse

from src.utils.general import read_subfamilies, create_directory, merge_commflags_with_kwargs, valid_file, parse_args
from src.utils.msa import load_msa
from src.grouping.upgma_subfam_grouping import run as run_upgma
from src.analysis.shannon_entropy import shannon_entropy_multiprocess,superfamily_shannon_entropy
from src.io.output import write_references
from src.io.input import read_temp_tea_file


ARGUEMENTS = {
	'm':{
		'flag':'msa_file',
		'required':True,
		'help':"File path to Multiple Sequence Alignment",
		'type':valid_file,
	},
	's': {
		'flag': 'subfamilies_file',
		'help': "File location of subfamily assignment (TXT file)",
		'type': valid_file,
	},
	'o': {
		'flag': 'outdir',
		'help': "Output file location",
		'type': str,
		'default': "TEA"
	},
	'p': {
		'flag': 'program_mode',
		'help': "Run program in standard TEA or objective TEAO",
		'type': str,
		'choices': ['TEA', 'TEAO'],
		'default': "TEAO"
	},
	'r': {
		'flag': 'reference_id',
		'help': "Reference ID(s)",
		'type': str,
	},
	'x': {
		'flag': 'similarity_matrix',
		'help': "Calculate scores according to provided matrix",
		'type': str,
	},
	't': {
		'flag': 'threads',
		'help': "Number of threads",
		'type': int,
		'default': 1
	}
}

def TEA(subfamilies_file:str,outdir:str,msa:pd.DataFrame,threads:int,mode:str,msa_file:str,file_prefix:str="TEA") -> np.ndarray:

	print("\n\tStarting subfamily entropy calculations...")

	temp_entropy_file_name = "intermediate_average_entropies.teao"
	if mode == "TEAO":
		subfamilies_file,_ = run_upgma(msa=msa,threads=threads,outdir=outdir)
		temp_entropy_file_name = "intermediate_average_entropies.teao"
	
	## Default to TEA-O messages
	if mode == "TEA" and not subfamilies_file:
		print("\n\t\t\t[W] Subfamily definition file not provided, defaulting to TEA-O algorithm...")
	if mode == "TEA" and not os.path.isfile(subfamilies_file):
		print(f"\n\t\t\t[W] Subfamily definition file {subfamilies_file} not found, defaulting to TEA-O algorithm...")
		subfamilies_file,_ = run_upgma(msa_file=msa_file,threads=threads,outdir=outdir)

	subfamilies:dict = read_subfamilies(subfamilies_file)

	temp_entropy_file = f"{outdir}/.data/{temp_entropy_file_name}"

	start = time.time()

	_,num_of_MSA_res = msa.loc[msa.index[0],'array'].shape

	keys = [x for x in subfamilies.keys()]

	results = {}

	if os.path.isfile(temp_entropy_file):
		print(f"\n\t\tPrevious subfamily entropies found, resuming at checkpoint...")
		results = read_temp_tea_file(temp_entropy_file)

	shannon_entropy_multiprocess(msa=msa,temp_entropy_file=temp_entropy_file,results=results,subfamilies=subfamilies,keys=keys,threads=threads)

	results = read_temp_tea_file(temp_entropy_file)

	E_i = np.zeros((1,num_of_MSA_res))
	for index in sorted(results.keys()):
		E_i += results[index]
	E_i /= len(results.keys())

	print(f"\n\t\tTEA-O processed {len(keys)} in {time.time()-start:.0f}s")

	with open(f"{outdir}/{file_prefix}.average_entropy",'w') as OUT:
		for msa_pos,se in enumerate(E_i[0]):
			OUT.write(f"{msa_pos}\t{se:.4f}\n")

	return E_i

def run(args:argparse.Namespace|None=None,**kwargs) -> None:

	"""
	Performs the (T)wo (E)ntropy (A)nalysis analysis on a given MSA, finding the global entropy,
	and the average entropy across subfamilies.

	Parameters:
	msa_file (str): File containing MSA for protein of interest
	subfamily_file (str): File containing user defined subfamily groupings for proteins in the provided MSA
	
	Returns:
	None
	"""

	args = merge_commflags_with_kwargs(
		cli_args=args,
		function_args=ARGUEMENTS,
		**kwargs
	)

	print(f"\n\tTwo Entropy Analysis started")
	
	create_directory(
		args.outdir
	)

	msa:pd.DataFrame = load_msa(
		args.msa_file,
		outdir=args.outdir
	)
	
	gE_i,a_i,G_i = superfamily_shannon_entropy(
		msa=msa,
		outdir=args.outdir
	)

	aE_i:np.array = TEA(
		subfamilies_file=args.subfamilies_file,
		outdir=args.outdir,
		msa=msa,
		threads=args.threads,
		mode=args.program_mode,
		msa_file=args.msa_file
	)

	write_references(
		reference_id=args.reference_id,
		msa=msa,
		global_SE=gE_i,
		average_SE=aE_i,
		residues=a_i,
		gaps=G_i,
		outdir=args.outdir
	)

	return None

if __name__ == '__main__':

	run(args=parse_args(function_args=ARGUEMENTS))