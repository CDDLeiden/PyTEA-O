import pathlib
import tempfile
import importlib
import argparse
from matplotlib import pyplot as plt

from PyTEAO import Tree, PhyloTree, UserTree
from PyTEAO import MSA
from PyTEAO import TwoEntropyAnalysis as TEA


def create_lineage_temp_file(nodes:list) -> pathlib.Path:

	temp_file,file_path = tempfile.mkstemp()
	temp_path = pathlib.Path(file_path)

	with open(temp_file,'w') as tmp:

		node:Tree.Node
		for id,node in enumerate(nodes):
			tmp.write(f"{id}\t1\t{",".join(node.accessions)}\n")

	return temp_path


def extract_entropy_values(tea:TEA):

	reference_indicies = tea.msa.get_reference_indices(tea.msa.reference_accession)

	global_entropy_values = tea.global_entropy[reference_indicies]
	global_entropy_values = (global_entropy_values - global_entropy_values.min()) / (global_entropy_values.max() - global_entropy_values.min())
	average_entropy_values = tea.average_entropy[reference_indicies]
	average_entropy_values = (average_entropy_values - average_entropy_values.min()) / (average_entropy_values.max() - average_entropy_values.min())
	
	return global_entropy_values, average_entropy_values


def nodes_at_num_subgroups(tree:Tree) -> dict:
	nodes_by_distance = tree.nodes_by_distance
	return {len(nodes_by_distance[dist]):nodes_by_distance[dist] for dist in nodes_by_distance.keys()}


def plot_validation_figure(tea:TEA, outdir:pathlib.Path, plot_title:str, plot_name:str) -> None:

	global_entropy_values, average_entropy_values = extract_entropy_values(tea)

	plt.figure(figsize=(4,4))
	plt.plot(average_entropy_values,global_entropy_values,marker='.',linestyle='None',markersize=6,color='blue',linewidth=0.05)
	plt.xlabel('Average Entropy')
	plt.ylabel('Global Entropy')
	plt.xticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
	plt.yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
	plt.title(plot_title)
	plt.margins(0)
	plt.subplots_adjust(left=0.15, bottom=0.15, right=0.98, top=0.95)
	plt.savefig(outdir/plot_name,dpi=300,bbox_inches='tight',pad_inches=0.05)


def validate(args:argparse.Namespace) -> None:

	outdir = pathlib.Path("./ValidationResults")

	validation_msa = importlib.resources.files("PyTEAO.package_data.validation_data")/"Ye-et-al_2008_alignment.fasta"

	msa = MSA(validation_msa,outdir=outdir)

	phylo_tree = PhyloTree(msa)

	nodes_by_num_subgroups = nodes_at_num_subgroups(phylo_tree)

	subgroups = [[25,'A'],[5,'B']]

	for datum in subgroups:

		groups = datum[0]

		nodes = [phylo_tree._nodes[node_id] for node_id in nodes_by_num_subgroups[groups]]

		temp_lineage_file = create_lineage_temp_file(nodes)

		tree = UserTree(msa,temp_lineage_file)

		plot_validation_figure(
			TEA(msa,tree,outdir=outdir/'TEA'),
			outdir/'TEA',
			f"Fig 2{datum[1]}: TEA-O Ye et al 2008 ({groups} subgroups)",
			f"Fig2{datum[1]}_TEA_plot_S{groups}.png"
		)

		temp_lineage_file.unlink()

	plot_validation_figure(
		TEA(msa,phylo_tree,outdir=outdir/'TEA'),
		outdir/'TEA',
		f"Fig 2C: TEA-O Ye et al 2008",
		f"Fig2C_TEA-P_plot.png"
	)