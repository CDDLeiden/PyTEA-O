#!/usr/bin/env python3

from argparse import ArgumentParser,Namespace
from src.analysis.tea import run as tea
from src.utils.general import get_file_name
from os.path import dirname,abspath
# from protein_descriptors import run as protein_descriptors
# from protein_descriptors import get_file_name
# from graph_results import run as graph_results

def run(args:Namespace) -> None:

	print("\n\n\tPerforming Two-Entropy Analysis.")
	tea(
		msa_file=args.msa_file,
		subfamilies_file=args.subfamily_file,
		outdir=args.outdir,
		mode=args.tea_mode,
		reference_id=args.reference,
		threads=args.threads
	)

	# print("\n\n\tIdentifying physiochemical patterns of residue positions.")
	# protein_descriptors(
	# 	outdir=args.outdir,
	# 	descriptor_file=args.descriptors,
	# 	summary_file=f"{args.outdir}/shannon_entropy.summary",
	# )

	# print("\n\n\tGraphing TEA-O results.")
	# graph_results(
	# 	shannon_entropy_summary_file=f"{args.outdir}/shannon_entropy.summary",
	# 	outdir=f"{args.outdir}",
	# 	highlight_residue_file=args.highlight_file,
	# 	subset_file=args.subset_file,
	# 	descriptor_file=f"{args.outdir}/prot_descriptors_{get_file_name(args.descriptors)}.tsv",
	# 	configuration_file=args.configuration
	# )


	print("\n\n\tTEA-O analysis has successfully completed.")
	return None

def parse_args() -> Namespace:

	parser = ArgumentParser()

	parser.add_argument("-m","--msa_file",required=True,type=str)
	parser.add_argument("-o","--outdir",default="TEA_ANALYSIS",type=str)
	parser.add_argument("-r","--reference",default=None,type=str)
	parser.add_argument("-e","--tea_mode",default="TEAO",choices=["TEA","TEAO"])
	parser.add_argument("-t","--threads",type=int,default=1)
	parser.add_argument("-d","--descriptors",default=f"{dirname(abspath(__file__))}/../DB_files/ZscaleSandberg.txt")
	parser.add_argument("-l","--highlight_file",default=None,type=str)
	parser.add_argument("-f","--subfamily_file",default=None,type=str)
	parser.add_argument("-s","--subset_file",default=None,type=str)
	parser.add_argument("-c","--configuration",default=f"{dirname(abspath(__file__))}/../config/graph_results.cfg")

	return parser.parse_known_args()[0]

if __name__ == "__main__":

	args = parse_args()
	run(args)