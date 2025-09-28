#!/usr/bin/env python3

import argparse

import src.utils.general as gen_util
import src.utils.msa as MSA
import src.analysis.phylotree as PhyloTree
import src.analysis.twoentropyanalysis as TEA
from src.visualization import plotmanager as PlotManager

# def parse_args() -> Namespace:

# 	parser = ArgumentParser()

# 	parser.add_argument("-m","--msa_file",required=True,type=str)
# 	parser.add_argument("-o","--outdir",default="TEA_ANALYSIS",type=str)
# 	parser.add_argument("-r","--reference",default=None,type=str)
# 	parser.add_argument("-p","--program_mode",default="TEAO",choices=["TEA","TEAO"])
# 	parser.add_argument("-t","--threads",type=int,default=1)
# 	parser.add_argument("-d","--descriptors",default=f"{dirname(abspath(__file__))}/DB_files/ZscaleSandberg.txt")
# 	parser.add_argument("-l","--highlight_file",default=None,type=str)
# 	parser.add_argument("-f","--subfamily_file",default=None,type=str)
# 	parser.add_argument("-s","--subset_file",default=None,type=str)
# 	parser.add_argument("-c","--configuration",default=f"{dirname(abspath(__file__))}/config/graph_results.cfg")

# 	return parser.parse_known_args()[0]

ARGUEMENTS = {
	'm':{
		'flag':'msa_file',
		'required':True,
		'help':"File path to Multiple Sequence Alignment",
		'type':gen_util.valid_file,
		'CLI': True,
		'call': False
	},
	'o': {
		'flag': 'outdir',
		'help': "Output file location",
		'type': str,
		'default': "TEA",
		'type': gen_util.valid_directory,
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
	'p': {
		'flag': 'program_mode',
		'help': "Run program in standard TEA or objective TEAO",
		'type': str,
		'choices': ['TEA', 'TEAO'],
		'default': "TEAO",
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
	},
	's': {
		'flag': 'subfamilies_file',
		'help': "File location of subfamily assignment (TXT file)",
		'type': gen_util.valid_file,
		'CLI': True,
		'call': True
	},
	'l': {
		'flag': 'plot_layout',
		'help': "Command line override for modifying subplots produced in figure",
		'type': str,
		'CLI': True,
		'call': True
	}
}


def run(args:argparse.Namespace|None=None,**kwargs) -> None:

	args:argparse.Namespace = gen_util.merge_commflags_with_kwargs(
		cli_args=args,
		function_args=ARGUEMENTS,
		**kwargs
	)

	print("Running MSA")
	msa = MSA.MSA(args.msa_file,args.outdir,threads=args.threads)

	print("Running Tree")
	tree = PhyloTree.PhyloTree(msa,threads=args.threads)

	print("Running TEA")
	tea = TEA.TwoEntropyAnalysis(msa,tree,threads=args.threads,outdir=args.outdir/"TEA")

	print("Plotting TEA")
	PlotManager.PlotManager(tea=tea,subplots='111111',outdir=args.outdir/"plots")


	return None

if __name__ == "__main__":

	run(args=gen_util.parse_args(function_args=ARGUEMENTS))