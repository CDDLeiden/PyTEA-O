#!/usr/bin/env python3

from argparse import ArgumentParser
from convert_alignment import run as convert_alignment
from upgma_subfam_grouping import run as upgma_subfam_grouping
from shannon_entropy import run as shannon_entropy
from z_scales import run as z_scales
from graph_results import run as graph_results

def run(args:ArgumentParser.parse_known_args) -> None:

	file_prefix = args.msa_file.split("/")[-1].split(".")[0]


	print("\n\n\tConverting MSA into compatible one-liner format.")
	convert_alignment(
		msa_file=args.msa_file,
		outdir=args.outdir,
		reference=args.reference,
		desired_alignment_type='one_liner',
		extract_fastas=True
	)

	print("\n\n\tGenerating UPGMA-based subfamily groupings.")
	upgma_subfam_grouping(
		msa_file=args.msa_file,
		outdir=args.outdir
	)

	print("\n\n\tPerforming Two-Entropy Analysis.")
	shannon_entropy(
		path_msa=args.msa_file,
		subfamilies_file=f"{args.outdir}/{file_prefix}.subfamilies",
		outdir=args.outdir,
		mode=args.tea_mode,
		reference_id=args.reference,
		threads=args.threads
	)

	exit()

	print("\n\n\tIdentifying physiochemical patterns of residue positions.")
	z_scales(
		msa_file=args.msa_file,
		outdir=args.outdir,
		descriptors=args.descriptors
	)

	print("\n\n\tGraphing TEA-O results.")
	graph_results(
		shannon_entropy_file=f"{args.outdir}/{args.reference}.shannon_entropy.txt",
		outdir=f"{args.outdir}",
		highlight_file=args.highlight_file,
		concensus_sequence_file=f"{args.outdir}/consensus_logo.txt",
		average_entropy_file=f"{args.outdir}/teao_average_shannon_entropy.txt",
		subset_file=args.subset_file,
		zscale_file=f"{args.outdir}/zscales.txt",
		configuration=args.configuration
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
	GetOptions.add_argument("-l","--highlight_file")
	GetOptions.add_argument("-s","--subset_file")
	GetOptions.add_argument("-c","--configuration",default=f"{dirname(abspath(__file__))}/../config/graph_results.cfg")

	run(GetOptions.parse_known_args()[0])