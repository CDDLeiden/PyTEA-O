#!/usr/bin/env python3

import argparse

import src.utils.general as gen_util
import src.utils.msa as MSA
import src.analysis.phylotree as PhyloTree
import src.analysis.twoentropyanalysis as TEA
from src.visualization import plotmanager as PlotManager

ARGUEMENTS = {
	'm':{
		'flag':'msa_file',
		'required':True,
		'help':"File path to Multiple Sequence Alignment",
		'type':gen_util.valid_file,
	},
	'o': {
		'flag': 'outdir',
		'help': "Output file location",
		'type': str,
		'default': "PyTEA-O_Results",
		'type': gen_util.valid_directory,
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
		'type': gen_util.valid_file,
	},
	'g': {
		'flag': 'highlight_file',
		'help': "File location for residues to highlight",
		'type': gen_util.valid_file,
	},
	'l': {
		'flag': 'plot_layout',
		'help': "Command line override for modifying subplots produced in figure",
		'type': gen_util.binary,
		'default': '111111'
	}
}

def run(args:argparse.Namespace|None=None) -> None:

	args:argparse.Namespace = gen_util.parse_args(
		function_args=ARGUEMENTS,
	)

	msa = MSA.MSA(args.msa_file,args.outdir,threads=args.threads)

	tree = PhyloTree.PhyloTree(msa,threads=args.threads)

	tea = TEA.TwoEntropyAnalysis(msa,tree,threads=args.threads,outdir=args.outdir/"TEA")

	figure = PlotManager.PlotManager(tea=tea,subplots=args.plot_layout,outdir=args.outdir/"plots",highlight_file=args.highlight_file)

	figure.save_fig(file_type="png")


if __name__ == "__main__":

	run(args=gen_util.parse_args(function_args=ARGUEMENTS))