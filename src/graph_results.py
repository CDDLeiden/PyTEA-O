#!/usr/bin/env python3

from matplotlib import pyplot as plt
from random import randint
import matplotlib.colors as mcolors
from matplotlib.ticker import MultipleLocator, AutoMinorLocator
from os.path import exists
from os import makedirs
from statistics import median, mean
import numpy as np
from math import log

def __read_consensus_file(consensus_file:str=None,data:dict=None,msa_links:dict=None) -> list:

	with open(consensus_file,'r') as IN:
		for line in IN:
			line = line.strip()

			if line == "":
				continue

			if line[0] == "#":
				continue

			msa_pos,res_nums = line.split()

			msa_pos = int(msa_pos)

			if msa_pos in msa_links.keys():
				non_gapped_data = [x for x in res_nums.split(";") if x.split(":")[0]]
				data[msa_links[msa_pos]]["CON_RES"] = non_gapped_data[0].split(":")[0]
				data[msa_links[msa_pos]]["RES_COUNTS"] = {x.split(":")[0]:float(x.split(":")[-1]) for x in non_gapped_data}
				data[msa_links[msa_pos]]["NUM_RES"] = len(non_gapped_data)
				temp_data = [int(x.split(":")[-1]) for x in non_gapped_data]
				data[msa_links[msa_pos]]["DEV"] = abs(mean(temp_data)-median(temp_data))/(mean(temp_data))

	return data

def __read_res_highlight_file(highlight_file:str=None) -> dict:

	highlights = {}
	if exists(highlight_file):
		source = None
		with open(highlight_file,'r') as IN:
			for line in IN:
				line = line.strip()

				if line == "":
					continue

				if line[0] == "#":
					source = line.replace("#","").strip()
					highlights[source] = []
					continue

				values = line.split("-")

				if len(values) == 2:
					for x in range(int(values[0]),int(values[1])+1):
						if int(x) not in highlights[source]:
							highlights[source].append(int(x))
				else:
					if int(values[0]) not in highlights[source]:
						highlights[source].append(int(values[0]))

	return highlights

def __read_SE_data(SE_file:str=None) -> list:

	if not exists(SE_file):
		print(f"\n\n\t[E]  {SE_file} does not exist!\n")
		exit

	data = {}
	msa_links = {}

	with open(SE_file,'r') as IN:
		for line in IN:
			line = line.strip()

			if line == "":
				continue

			if line[0] == "#":
				continue

			pos,msa_pos,res,se,fse,num_seqs,fog = line.split()

			data[int(pos)+1] = {"RES":res,"SE":abs(float(fse)),"NUM_RES":None,"MSA_POS":msa_pos}
			msa_links[int(msa_pos)] = int(pos) + 1

	return data,msa_links


def plot(SE_file:str = None, consensus_file:str = None, output_dir: str = 'SE_graphics', highlight_file:str = None) -> None:

	aas = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']

	# Read in Shannon Entropy
	data,msa_links = __read_SE_data(SE_file=SE_file)

	# Read in Concensus Sequence
	data = __read_consensus_file(consensus_file=consensus_file,data=data,msa_links=msa_links)
	
	## Read in residues to highlight in the graphic
	highlights = __read_res_highlight_file(highlight_file=highlight_file)
	
	# Reference residue sequences
	res_nums = sorted(data.keys())

	# Number of sequences present in MSA
	num_of_seq = sum(x for x in data[res_nums[0]]["RES_COUNTS"].values())

	# Number of residues in the reference sequence
	len_of_seq = len(res_nums)
	# Font size of tick labels
	font_size = 6
	# Width of the graphic as to not have labels overlapping each other
	width = font_size*(1/72)*(len_of_seq)*1.5
	
	####################################################################################################################
	## Setting up graphics
	####################################################################################################################

	fig = plt.figure()
	
	## Specifying plot placements
	# Shannon Entropy 1st
	se_graph = plt.subplot2grid((5,1),(0,0))
	# Residue presense matrix 2nd, spans 2 rows
	res_mat_graph = plt.subplot2grid((5,1),(1,0),rowspan=2)
	# Number of different residues present at site
	num_ress = np.zeros(len_of_seq)
	res_num_graph = plt.subplot2grid((5,1),(3,0))
	# Number of sequences belonging to that residue
	res_dist_graph = plt.subplot2grid((5,1),(4,0))
	
	## Load CSS colors for iterative coloring
	colors = mcolors.XKCD_COLORS
	color_keys = [
		"aqua","blue","chartreuse","coral","crimson",
		"darkgreen","fuchsia","goldenrod","indigo","navy",
		"olive","orange","orangered","orchid","salmon",
		"teal","tomato","yellowgreen","grey","khaki"
	]

	####################################################################################################################
	## Shannon Entropy Graphing
	####################################################################################################################

	# Shannon Entropies for the reference
	se_values = [data[x]["SE"] for x in res_nums]

	# Residue and residue number labels
	labels = [f"{data[x]["RES"]} [{x}]"  if x%2 != 0 else f"{data[x]["RES"]}" for x in res_nums]

	## Expand the graph slightly past the plotted data
	se_graph.set_xlim([(res_nums[0]-0.5),res_nums[-1]+0.5])
	se_graph.set_ylim([-0.05,1.05])

	# Add separator lines to make residue visualization easier
	for x in res_nums:
		se_graph.plot([x-0.5,x-0.5],[-0.5,1.05],color='w',linestyle='-',linewidth=0.75)
	
	## Add "Highlighted Entropy Zones"
	# Green (0.0-0.25)
	# Yellow (0.25-0.75)
	# Red (0.75-1.0)
	se_graph.add_patch(plt.Rectangle((res_nums[0]-1,-0.05),(res_nums[-1]+2),0.3,color='g',alpha=0.1))
	se_graph.add_patch(plt.Rectangle((res_nums[0]-1,.25),(res_nums[-1]+2),0.5,color='y',alpha=0.1))
	se_graph.add_patch(plt.Rectangle((res_nums[0]-1,.75),(res_nums[-1]+2),1,color='r',alpha=0.1))
	
	# Highlight residues with near equal representation in MSA
	for i,val_x in enumerate(res_nums):
		num_res = data[val_x]["NUM_RES"]
		num_res -= 1 if '-' in data[val_x]["RES_COUNTS"].keys() else 0
		num_res -= 1 if 'X' in data[val_x]["RES_COUNTS"].keys() else 0
		num_ress[i] = num_res
		if num_res < 10:
			continue
		se_graph.plot(val_x,data[val_x]["SE"],marker="*",markersize=10*(1-data[val_x]["DEV"]),color='b')

	## Highlight user-given residues
	if exists(highlight_file):
		used_colors = []
		for source in highlights.keys():
			color = colors[f"xkcd:{color_keys[randint(0,len(color_keys)-1)]}"]
			while color in used_colors:
				if sorted(used_colors) == sorted(color_keys):
					used_colors = []
			for val_x in highlights[source]:
				se_graph.add_patch(plt.Rectangle((val_x-0.5,-0.05),1,1.2,facecolor=color,edgecolor="None"))
				used_colors.append(color)

	## Plot the data
	se_graph.plot(res_nums,se_values,linewidth=0.5,color='k')

	## Prettify the Y-axis
	se_graph.yaxis.set_ticks([i*0.2 for i in range(int(1/0.2)+1)])
	se_graph.yaxis.set_ticklabels([0.0,0.2,0.4,0.6,0.8,1.0])

	se_graph.tick_params(axis='x',which='both',bottom=False,labelbottom=False)
	se_graph.xaxis.grid(color='w',linestyle='-',linewidth=0.75,which='minor')
	se_graph.spines['bottom'].set_visible(False)
	axt = se_graph.secondary_xaxis('top')
	axt.set_xticks(res_nums)
	axt.set_xticklabels(labels,font='monospace',fontsize=font_size,ha='center',va='bottom',rotation='vertical')

	####################################################################################################################
	## Residue Presence Matrix
	####################################################################################################################

	## Get a positive/negative matrix of what residues are present in the MSA at each residue
	present_res = np.zeros((len(aas),len(res_nums)))
	for row_index,aa in enumerate(aas):
		for col_index,val_x in enumerate(res_nums):
			if aa in data[val_x]["RES_COUNTS"].keys():
				present_res[row_index][col_index] = 1 if data[val_x]["RES_COUNTS"][aa] > 0 else 0

	res_mat_graph.matshow(present_res,aspect='auto',cmap='tab20_r')

	## Setup Y-Axis
	res_mat_graph.set_yticks([i for i in range(len(aas))],labels=aas,font='monospace')
	res_mat_graph.yaxis.set_minor_locator(MultipleLocator(1,offset=-0.5))
	res_mat_graph.yaxis.grid(color='w',linestyle='-',linewidth=0.75,which='minor')
	res_mat_graph.tick_params(axis='y',which='minor',left=False)

	## Setup X-Axis
	# Set major tick for each residue [0 to pos-1]
	res_mat_graph.xaxis.set_ticks([i-1 for i in res_nums])
	# Set and rotate major tick labels
	labels = [f"{val_x: >{len(str(len_of_seq))+1}}" if val_x%2 != 0 else "" for val_x in res_nums]
	res_mat_graph.set_xticklabels(labels=labels,rotation='vertical',font='monospace',fontsize=font_size,va='center',ha='center')
	# Show major tick labels, but not the ticks themselves, and only on the top
	res_mat_graph.tick_params(axis='x',which='major',top=False,bottom=False,labeltop=True,labelbottom=False)

	# Set minor ticks for each residue, but offset by -0.5 to create fake gridlines
	res_mat_graph.xaxis.set_minor_locator(MultipleLocator(1,offset=-0.5))
	# Don't show any information regarding the minor ticks
	res_mat_graph.tick_params(axis='x',which='minor',top=False,bottom=False,labeltop=False,labelbottom=False)
	
	# Hide the top and bottom axes
	res_mat_graph.spines['top'].set_visible(False)
	res_mat_graph.spines['bottom'].set_visible(False)

	# Turn on X-axis gridlines using the minor ticks
	res_mat_graph.xaxis.grid(color='w',linestyle='-',linewidth=0.75,which='minor')

	####################################################################################################################
	## Number of Residue Bar Graph
	####################################################################################################################

	res_num_graph.bar(res_nums,num_ress,1)

	## Setup X-Axis
	# Set major tick for each residue [0 to pos-1]
	res_num_graph.xaxis.set_ticks([i for i in res_nums])
	# Set and rotate major tick labels
	res_num_graph.set_xticklabels(labels=labels,rotation='vertical',font='monospace',fontsize=font_size,va='center',ha='center')
	res_num_graph.tick_params(axis='x',which='major',top=False,bottom=False,labelbottom=False,labeltop=True)

	res_num_graph.xaxis.set_minor_locator(MultipleLocator(1,offset=-0.5))
	res_num_graph.grid(color='gainsboro',linestyle='-',linewidth=0.75,which='minor')
	res_num_graph.tick_params(axis='x',which='minor',top=False,bottom=False)

	res_num_graph.set_xlim([(res_nums[0]-0.5),res_nums[-1]+0.5])

	## Setup Y-Axis
	res_num_graph.set_ylim([-0.05,21])
	res_num_graph.set_yticks([0,5,10,15,20])
	
	# Show horizontal guidelines
	res_num_graph.plot([res_nums[0]-0.5,res_nums[-1]+0.5],[5,5],color='w',linestyle='--',linewidth=1)
	res_num_graph.plot([res_nums[0]-0.5,res_nums[-1]+0.5],[10,10],color='aqua',linestyle='-',linewidth=1)
	res_num_graph.plot([res_nums[0]-0.5,res_nums[-1]+0.5],[15,15],color='w',linestyle='--',linewidth=1)
	res_num_graph.plot([res_nums[0]-0.5,res_nums[-1]+0.5],[20,20],color='aqua',linestyle='-',linewidth=1)

	# Hide the top and bottom axes
	res_num_graph.spines['top'].set_visible(False)
	res_num_graph.spines['bottom'].set_visible(False)
	

	# ###################################################################################################################
	# ## Residue Distribution Bar Graph
	# ###################################################################################################################

	# Stacking bars require a y-value to place the bottom of the bar
	
	bottoms = np.zeros(len_of_seq)

	for i in range(len(aas)):
		res_vals = np.zeros(len_of_seq)
		for index,val_x in enumerate(res_nums):
			if i < len(data[val_x]["RES_COUNTS"].keys()):
				sorted_res_counts = sorted(data[val_x]["RES_COUNTS"].keys(),reverse=True,key=lambda x: data[val_x]["RES_COUNTS"][x])
				res_vals[index] = data[val_x]["RES_COUNTS"][sorted_res_counts[i]]/num_of_seq if sorted_res_counts[i] != "-" else 0
			else:
				res_vals[index] = 0
		res_dist_graph.bar(res_nums,res_vals,1,bottom=bottoms)
		bottoms += res_vals

	## Setup top X-Axis
	res_dist_graph.set_xlim([(res_nums[0]-0.5),res_nums[-1]+0.5])
	# Set major tick for each residue [0 to pos-1]
	res_dist_graph.xaxis.set_ticks([i for i in res_nums])
	# Set and rotate major tick labels
	res_dist_graph.set_xticklabels(labels=labels,rotation='vertical',fontsize=font_size,va='center',ha='center')
	res_dist_graph.tick_params(axis='x',which='major',top=False,bottom=False,labelbottom=False,labeltop=True)
	res_dist_graph.tick_params(axis='x',which='minor',top=False,bottom=False,labeltop=False,labelbottom=False)
	res_dist_graph.xaxis.set_minor_locator(MultipleLocator(1,offset=-0.5))
	res_dist_graph.grid(color='gainsboro',linestyle='-',linewidth=0.75,which='minor')
	## Setup bottom X-Axis
	rdg_bt = res_dist_graph.secondary_xaxis("bottom")
	rdg_bt.set_xticks([i for i in res_nums])
	rdg_bt.tick_params(which='both',bottom=False)
	labels = [f"[{x}] {data[x]["CON_RES"]}" if x%2 != 0 else data[x]["CON_RES"] for x in res_nums]
	rdg_bt.set_xticklabels(labels=labels,rotation='vertical',font='monospace',fontsize=font_size)

	## Setup Y-Axis
	res_dist_graph.yaxis.set_ticks([i*0.2 for i in range(int(1/0.2)+1)])
	res_dist_graph.yaxis.set_ticklabels([0.0,0.2,0.4,0.6,0.8,1.0])

	res_dist_graph.spines['top'].set_visible(False)
	
	## Saving Files

	file_prefix = SE_file.split("/")[-1].split(".")[0]
	makedirs(output_dir,mode=0o755,exist_ok=True)
	fig.set_size_inches(width,10)
	plt.savefig(f"{output_dir}/{file_prefix}.png",dpi=500)
	
	return None

def run(args=None) -> None:

	plot(SE_file=args.shannon_entropy_file,consensus_file=args.consensus_sequence_file,highlight_file=args.highlight_residues,output_dir=args.outdir)

	return None

if __name__ == "__main__":

	from argparse import ArgumentParser

	GetOptions = ArgumentParser()

	GetOptions.add_argument("-s","--shannon_entropy_file",required=True,type=str)
	GetOptions.add_argument("-o","--outdir",required=False,type=str,default="se_graphics")
	GetOptions.add_argument("-y","--highlight_residues",required=False,type=str)
	GetOptions.add_argument("-c","--consensus_sequence_file",required=False,type=str)
	GetOptions.add_argument("-a","--average_entropy_file",required=False)

	run(GetOptions.parse_known_args()[0])