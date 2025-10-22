#!/usr/bin/env python3

import os
import pathlib
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from PyTEAO import PhyloTree, TaxonTree, Tree, UserTree
from PyTEAO import TwoEntropyAnalysis as TEA
from PyTEAO import MSA

def subgroup_by_upgma_split(tree:Tree,nodes_by_distance:dict,num_subgroups:int, output_dir:pathlib.Path) -> pathlib.Path|None:

	outfile = output_dir/f"subgroups_upgma_split_{num_subgroups}.headers"

	for distance in sorted(nodes_by_distance.keys()):

		if len(nodes_by_distance[distance]) == num_subgroups:

			with outfile.open('w') as OUT:

				header = ['##', 'Node_id', 'Branch_distance', 'Accessions']
				OUT.write(f"{'\t'.join(header)}\n")
				count = 0
				for node in nodes_by_distance[distance]:
					accessions = tree._nodes[node].accessions
					OUT.write(f"{count}\t{count}\t{','.join(accessions)}\n")
					count += 1

			return outfile
			
	return None


def run(msa_file:pathlib.Path, outdir:pathlib.Path, threads:int, plot_title, plot_name,subfamilies_file:pathlib.Path=None) -> None:

	## Load the MSA
	msa = MSA(msa_file,outdir,threads=threads)

	## Build Tree for Two Entropy Calculations
	tree:Tree
	if subfamilies_file is not None:
			tree = UserTree(msa,subfamilies_file)
	else:
		tree = PhyloTree(msa)

	## Perform Two Entropy Analaysis
	tea = TEA(msa,tree,threads=threads,outdir=outdir)

	reference_indicies = msa.get_reference_indices(msa.reference_accession)

	global_entropy_values = tea.global_entropy[reference_indicies]
	global_entropy_values = (global_entropy_values - global_entropy_values.min()) / (global_entropy_values.max() - global_entropy_values.min())
	average_entropy_values = tea.average_entropy[reference_indicies]
	average_entropy_values = (average_entropy_values - average_entropy_values.min()) / (average_entropy_values.max() - average_entropy_values.min())

	plt.figure(figsize=(4,4))
	plt.plot(average_entropy_values,global_entropy_values,marker='.',linestyle='None',markersize=6,color='blue',linewidth=0.05)
	plt.xlabel('Average Entropy')
	plt.ylabel('Global Entropy')
	plt.xticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
	plt.yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
	plt.title(plot_title)
	plt.margins(0)
	plt.subplots_adjust(left=0.15, bottom=0.15, right=0.98, top=0.95)
	plt.savefig(output_dir/plot_name,dpi=300,bbox_inches='tight',pad_inches=0.05)


if __name__ == "__main__":

	validation_directory = pathlib.Path(__file__).resolve(strict=False).parent

	msa_file = validation_directory/'Ye-et-al_2008_alignment.fasta'
	output_dir = validation_directory/'output'
	threads = 2

	from PyTEAO import MSA

	## Load the MSA
	msa = MSA(msa_file,output_dir,threads=threads)

	from PyTEAO import PhyloTree, TaxonTree, Tree

	## Build Tree for Two Entropy Calculations
	tree = PhyloTree(msa)

	nodes_by_distance = tree.nodes_by_distance

	# For TEA mode with UPGMA splits
	subgroup_by_upgma_split(tree,nodes_by_distance, 25, output_dir)
	subgroupfile_S25 = output_dir/'subgroups_upgma_split_25.headers'
	run(
		msa_file,
		outdir=output_dir,
		threads=threads, 
		subfamilies_file=subgroupfile_S25, 
		plot_title='Fig 2A: TEA-O Ye et al 2008 (25 subgroups)', 
		plot_name='Fig2A_TEA_plot_S25.png',
	)
	
	subgroup_by_upgma_split(tree, nodes_by_distance, 5, output_dir)
	subgroupfile_S5 = output_dir/'subgroups_upgma_split_5.headers'
	run(
		msa_file, outdir=output_dir,
		threads=threads,
		subfamilies_file=subgroupfile_S5,
		plot_title='Fig 2B: TEA-O Ye et al 2008 (5 subgroups)',
		plot_name='Fig2B_TEA_plot_S5.png',
	)

	# For TEA-O mode without subgroups
	run(
		msa_file,
		outdir=output_dir,
		threads=threads,
		plot_title='Fig 2C: TEA-O Ye et al 2008',
		plot_name='Fig2C_TEA-O_plot.png',
	)
	
	# Remove subgroup files
	subgroupfile_S25.unlink()
	subgroupfile_S5.unlink()