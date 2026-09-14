#!/usr/bin/env python3

import click
import typing
import pathlib

from PyTEAO import general, validation

@click.command()
@click.option(
	"-m","--msa-file",
	help="File containing aligned sequences",
	type=click.Path(path_type=pathlib.Path,dir_okay=False,readable=True,exists=True),
	required=True
)
@click.option(
	"-o","--outdir",
	help="Output directory to place results",
	type=click.Path(path_type=pathlib.Path,file_okay=False,readable=True,writable=True,exists=False),
	default="PyTEA-O",
	required=False
)
@click.option(
	"-r","--reference-id",
	help="Accession number of sequence in MSA to use as base in figure generation",
	type=str,
	default=None
)
@click.option(
	"-t","--threads",
	help="Number of threads allowed for process use",
	type=int,
	default=1
)
@click.option(
	"-s","--subfamilies-file",
	help="File containing subfamily assignments",
	type=click.Path(path_type=pathlib.Path,dir_okay=False,exists=True,readable=True),
	required=False
)
@click.option(
	"-g","--highlight-file",
	help="File containing residue-specific highlights",
	type=click.Path(path_type=pathlib.Path,dir_okay=False,exists=True,readable=True),
	required=False
)
@click.option(
	"-l","--plot-layout",
	help="Command line override for modifying subplots present in output figure",
	type=general.binary,
	default='111111',
	required=False
)
@click.option(
	"-b","--tree-type",
	help="Method for grouping sequences in the provided MSA",
	type=click.Choice(['dist','lin']),
	default='dist',
	required=False
)
@click.option(
	"-v","--validate",
	help="Run TEA-O validation analysis",
	type=bool,
	is_flag=True,
	default=False,
	required=False
)
def main(
		msa_file:pathlib.Path,
		*,
		outdir:pathlib.Path=pathlib.Path("./PyTEA-O"),
		reference_id:str|None=None,
		threads:int=1,
		subfamilies_file:pathlib.Path|None=None,
		highlight_file:pathlib.Path|None=None,
		plot_layout:general.binary='111111',
		tree_type:typing.Literal['dist','lin'] = 'dist',
		validate:bool=False
	) -> None:

	(outdir := pathlib.Path(outdir)).mkdir(mode=0o755,parents=True,exist_ok=True)

	if not validate and msa_file is None:
		raise ValueError(f"The --msa_file argument is required unless --validate is used.")


	## Validate pipeline against Ye's Results
	if validate:
		validation.validate(outdir)
		return
	

	from PyTEAO import MSA

	## Load the MSA
	msa = MSA(
			msa_file,
			outdir,
			threads=threads,
			reference_accession=reference_id
		)


	from PyTEAO import PhyloTree, TaxonTree, UserTree, Tree

	## Build Tree for Two Entropy Calculations
	tree:Tree
	if subfamilies_file is not None:
		tree = UserTree(msa,subfamilies_file)
	elif tree_type == 'lin':
		tree = TaxonTree(msa)
	elif tree_type == 'dist':
		tree = PhyloTree(msa)
	else:
		raise ValueError(f"Invalid tree_type {tree_type} provided, 'lineage' or 'distance_matrix' supported.")

	from PyTEAO import TwoEntropyAnalysis as TEA

	## Perform Two Entropy Analaysis
	tea = TEA(
			msa,
			tree,
			threads=threads,
			outdir=outdir,
		)

	from PyTEAO import PlotManager

	## Plot data
	figure = PlotManager(
			tea=tea,
			subplots=plot_layout,
			outdir=outdir,
			highlight_file=highlight_file
		)

	figure.save_fig(file_type="png")
	figure.save_fig(file_type="svg")


if __name__ == "__main__":

	main()