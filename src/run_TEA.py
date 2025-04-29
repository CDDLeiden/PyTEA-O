#!/usr/bin/env python3

from argparse import ArgumentParser
from shannon_entropy import run as shannon_entropy
from z_scales import run as z_scales
from graph_results import run as graph_results

def run(args:ArgumentParser.parse_known_args) -> None:

	file_prefix = args.msa_file.split("/")[-1].split(".")[0]

	print("\n\n\tPerforming Two-Entropy Analysis.")
	shannon_entropy(
		msa_file=args.msa_file,
		subfamilies_file=args.subfamily_file,
		outdir=args.outdir,
		mode=args.tea_mode,
		reference_id=args.reference,
		threads=args.threads
	)

	print("\n\n\tIdentifying physiochemical patterns of residue positions.")
	z_scales(
		msa_file=args.msa_file,
		outdir=args.outdir,
		descriptors=args.descriptors
	)

	print("\n\n\tGraphing TEA-O results.")
	graph_results(
		shannon_entropy_summary_file=f"{args.outdir}/shannon_entropy.summary",
		outdir=f"{args.outdir}",
		highlight_residue_file=args.highlight_file,
		subset_file=args.subset_file,
		zscale_file=f"{args.outdir}/zscales.tsv",
		configuration_file=args.configuration
	)


	print("\n\n\tTEA-O analysis has successfully completed.")
	return None

if __name__ == "__main__":

	from os.path import dirname,abspath

	GetOptions = ArgumentParser()

	GetOptions.add_argument("-m","--msa_file",required=True,type=str)
	GetOptions.add_argument("-o","--outdir",default="TEA_ANALYSIS",type=str)
	GetOptions.add_argument("-r","--reference",default=None,type=str)
	GetOptions.add_argument("-e","--tea_mode",default="TEAO",choices=["TEA","TEAO"])
	GetOptions.add_argument("-t","--threads",type=int,default=1)
	GetOptions.add_argument("-d","--descriptors",default='Zscale Sandberg')
	GetOptions.add_argument("-l","--highlight_file",default=None,type=str)
	GetOptions.add_argument("-f","--subfamily_file",default=None,type=str)
	GetOptions.add_argument("-s","--subset_file",default=None,type=str)
	GetOptions.add_argument("-c","--configuration",default=f"{dirname(abspath(__file__))}/../config/graph_results.cfg")

	run(GetOptions.parse_known_args()[0])