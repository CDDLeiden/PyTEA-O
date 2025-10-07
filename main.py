#!/usr/bin/env python3

import argparse

from PyTEAO import PhyloTree, TaxonTree, Tree, PlotManager, MSA
from PyTEAO import general,argparse
from PyTEAO import TwoEntropyAnalysis as TEA

ARGUEMENTS = {
	'm':{
		'flag':'msa_file',
		'required':True,
		'help':"File path to Multiple Sequence Alignment",
		'type':general.valid_file,
	},
	'o': {
		'flag': 'outdir',
		'help': "Output file location",
		'type': str,
		'default': "PyTEA-O_Results",
		'type': general.valid_directory,
	},
	'r': {
		'flag': 'reference_id',
		'help': "Reference ID",
		'type': str,
	},
	'p': {
		'flag': 'program_mode',
		'help': "Run program in standard TEA or objective TEAO",
		'type': str,
		'choices': ['TEA', 'TEAO'],
		'default': "TEAO",
	},
	't': {
		'flag': 'threads',
		'help': "Number of threads",
		'type': int,
		'default': 1,
	},
	's': {
		'flag': 'subfamilies_file',
		'help': "File location of subfamily assignment (TXT file)",
		'type': general.valid_file,
	},
	'g': {
		'flag': 'highlight_file',
		'help': "File location for residues to highlight",
		'type': general.valid_file,
	},
	'l': {
		'flag': 'plot_layout',
		'help': "Command line override for modifying subplots produced in figure",
		'type': general.binary,
		'default': '111111'
	},
	'b': {
		'flag': 'tree_type',
		'help': "Method to group sequences provided in the MSA",
		'type': str,
		'default': 'distance_matrix',
		'choices':['distance_matrix','lineage']
	}
}

def run(args:argparse.Namespace|None=None) -> None:

	## Create and parse arguments
	args:argparse.Namespace = general.parse_args(
		function_args=ARGUEMENTS,
	)

	## Load the MSA
	msa = MSA.MSA(args.msa_file,args.outdir,threads=args.threads)

	## Build Tree for Two Entropy Calculations
	tree:Tree
	if args.tree_type == 'lineage':
		tree = TaxonTree(msa)

	elif args.tree_type == 'distance_matrix':
		tree = PhyloTree(msa)

	else:
		raise ValueError(f"Invalid tree_type {args.tree_type} provided, 'lineage' or 'distance_matrix' supported.")


	## Perform Two Entropy Analaysis
	tea = TEA(msa,tree,threads=args.threads,outdir=args.outdir/"TEA")

	## Plot data
	figure = PlotManager(tea=tea,subplots=args.plot_layout,outdir=args.outdir/"plots",highlight_file=args.highlight_file)

	figure.save_fig(file_type="png")


if __name__ == "__main__":

	run(args=general.parse_args(function_args=ARGUEMENTS))