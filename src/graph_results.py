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

def plot_SE(SE_file:str = None, consensus_file:str = None, output_dir: str = 'SE_graphics', highlight_file:str = None) -> None:

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

			data[int(pos)] = {"RES":res,"SE":abs(float(fse)),"NUM_RES":None,"MSA_POS":msa_pos,"MED_RES":None}
			msa_links[int(msa_pos)] = int(pos)

	num_seqs = 0
	aas = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
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
				non_gapped_data = [x for x in res_nums.split(";") if x.split(":")[0] in aas]
				data[msa_links[msa_pos]]["RES_COUNTS"] = {x.split(":")[0]:float(x.split(":")[-1]) for x in non_gapped_data}
				num_seqs = sum([float(x.split(":")[-1]) for x in res_nums.split(";")])
				data[msa_links[msa_pos]]["NUM_RES"] = len(non_gapped_data)
				temp_data = [int(x.split(":")[-1]) for x in non_gapped_data]
				data[msa_links[msa_pos]]["DEV"] = abs(mean(temp_data)-median(temp_data))/(mean(temp_data))

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

	res_nums = sorted(data.keys())

	num_of_res = len(res_nums)
	font_size = 6
	width = font_size*(1/72)*(num_of_res)*1.5

	y = [data[x]["SE"] for x in res_nums]
	labels = [f"{data[x]["RES"]} - {x+1}"  if x%2 == 0 else f"{data[x]["RES"]}" for x in res_nums]

	fig = plt.figure()

	colors = mcolors.CSS4_COLORS
	color_keys = [x for x in colors]

	res_num_graph = plt.subplot2grid((5,1),(4,0))
	res_mat_graph = plt.subplot2grid((5,1),(1,0),rowspan=2)
	res_dist_graph = plt.subplot2grid((5,1),(3,0))

	############################################################
	## Shannon Entropy Graphing
	############################################################

	## Positional the plot in the subplot
	se_graph = plt.subplot2grid((5,1),(0,0))

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
	
	# # Highlight residues with near equal representation in MSA
	# for i,val_x in enumerate(x):
	# 	num_res = data[val_x]["NUM_RES"]
	# 	num_ress[i] = num_res
	# 	med_res = data[val_x]["MED_RES"]
	# 	if num_res < 10:
	# 		continue
	# se_graph.plot(val_x,data[val_x]["SE"],marker="*",markersize=10*(1-data[val_x]["DEV"]),color='b')

	## Highlight user-given residues
	if exists(highlight_file):
		used_colors = []
		for source in highlights.keys():
			color = colors[color_keys[randint(0,len(color_keys))]]
			while color in used_colors:
				if sorted(used_colors) == sorted(color_keys):
					used_colors = []
			for val_x in highlights[source]:
				se_graph.add_patch(plt.Rectangle((val_x-0.5,-0.05),1,1.2,facecolor=color,edgecolor="None"))
				used_colors.append(color)

	## Plot the data
	se_graph.plot(res_nums,y,linewidth=0.5,color='k')

	## Prettify the Y-axis
	se_graph.yaxis.set_ticks([i*0.2 for i in range(int(1/0.2)+1)])
	se_graph.yaxis.set_ticklabels(["",0.2,0.4,0.6,0.8,1.0])


	fig.set_size_inches(width,10)
	plt.savefig("MyTest.png",dpi=500)

	exit()

	# se_graph.tick_params(axis='x',which='both',bottom=False,labelbottom=False)
	# se_graph.xaxis.grid(color='w',linestyle='-',linewidth=0.75,which='minor')
	# se_graph.spines['bottom'].set_visible(False)
	# axt = se_graph.secondary_xaxis('top')
	# axt.set_xticks(x)
	# axt.set_xticklabels(labels,font='monospace',fontsize=font_size,ha='center',rotation='vertical')


	# ############################################################
	# ## Residue Distribution Bar Graph
	# ############################################################



	# num_ress = np.zeros(len(x))
	# res_dist_graph.bar(x,num_ress,1)
	# res_dist_graph.spines['top'].set_visible(False)
	# res_dist_graph.spines['bottom'].set_visible(False)
	# res_dist_graph.tick_params(axis='x',which='both',top=False,bottom=False,labeltop=False,labelbottom=False)
	# res_dist_graph.set_xlim([(x[0]-0.5),x[-1]+0.5])
	# res_dist_graph.xaxis.set_minor_locator(MultipleLocator(1,offset=-0.5))
	# res_dist_graph.yaxis.set_ticks([i*5 for i in range(int(20/5)+1)])
	# res_dist_graph.set_yticklabels(["",5,10,15,20])
	# res_dist_graph.grid(color='gainsboro',linestyle='-',linewidth=0.75,which='minor')
	# res_dist_graph.tick_params(axis='x',which='minor',top=False,bottom=False)
	# res_dist_graph.plot([x[0]-0.5,x[-1]+0.5],[10,10],color='aqua',linestyle='-',linewidth=2)
	# res_dist_graph.plot([x[0]-0.5,x[-1]+0.5],[20,20],color='aqua',linestyle='-',linewidth=2)
	# res_dist_graph.set_ylim([-0.05,21])
	
	# bottoms = np.zeros(len(x))
	# for i in range(len(aas)):
	# 	res_vals = np.zeros(len(x))
	# 	for index,val_x in enumerate(x):
	# 		if i < len(data[val_x]["RES_COUNTS"].keys()):
	# 			res_vals[index] = sorted(data[val_x]["RES_COUNTS"].values(),reverse=True)[i]/num_seqs if sorted(data[val_x]["RES_COUNTS"].keys(),reverse=True,key=lambda x: data[val_x]["RES_COUNTS"][x])[i] != "-" else 0
	# 		else:
	# 			res_vals[index] = 0
	# 	res_num_graph.bar(x,res_vals,1,bottom=bottoms)
	# 	bottoms += res_vals

	# ############################################################
	# ## Number of Residue Bar Graph
	# ############################################################

	# res_num_graph.xaxis.set_minor_locator(MultipleLocator(1,offset=-0.5))
	# res_num_graph.yaxis.set_ticks([i*0.2 for i in range(int(1/0.2)+1)])
	# res_num_graph.set_yticklabels(["",0.2,0.4,0.6,0.8,1.0])
	# res_num_graph.set_xlim([(x[0]-0.5),x[-1]+0.5])
	# res_num_graph.grid(color='gainsboro',linestyle='-',linewidth=0.75,which='minor')
	# res_num_graph.spines['top'].set_visible(False)
	# res_num_graph.xaxis.set_major_locator(MultipleLocator(25))
	# res_num_graph.tick_params(axis='x',which='major',top=False,labelbottom=True,labeltop=False,length=5)
	# res_num_graph.tick_params(axis='x',which='minor',top=False,bottom=False)
	# res_num_graph.set_ylim([0,1.05])

	# values = np.zeros((len(aas),len(x)))
	# for row_index,aa in enumerate(aas):
	# 	for col_index,val_x in enumerate(x):
	# 		if aa in data[val_x]["RES_COUNTS"].keys():
	# 			values[row_index][col_index] = 1 if data[val_x]["RES_COUNTS"][aa] > 0 else 0

	# ############################################################
	# ## Residue Presence Matrix
	# ############################################################

	# res_mat_graph.matshow(values,aspect='auto',cmap='tab20_r')
	# res_mat_graph.set_yticks([i for i in range(len(aas))],labels=aas,font='monospace')
	# res_mat_graph.yaxis.set_minor_locator(MultipleLocator(1,offset=-0.5))
	# res_mat_graph.yaxis.grid(color='w',linestyle='-',linewidth=0.75,which='minor')
	# res_mat_graph.xaxis.grid(color='w',linestyle='-',linewidth=0.75,which='minor')
	# res_mat_graph.xaxis.set_major_locator(MultipleLocator(25))
	# res_mat_graph.tick_params(axis='x',which='both',top=False,bottom=False,labeltop=False,labelbottom=False)
	# res_mat_graph.tick_params(axis='x',which='minor',labeltop=True,labelbottom=True)
	# res_mat_graph.xaxis.set_minor_locator(MultipleLocator(1,offset=-0.5))
	# res_mat_graph.set_xticklabels(minor=True,labels=[f" {val_x} " for val_x in x],rotation='vertical',fontsize=font_size,va='center',ha='center')
	# res_mat_graph.tick_params(axis='y',which='minor',left=False)
	# res_mat_graph.spines['top'].set_visible(False)
	# res_mat_graph.spines['bottom'].set_visible(False)

	# fig.subplots_adjust(hspace=0.15)
	# # plt.tight_layout()
	# file_prefix = SE_file.split("/")[-1].split(".")[0]
	# makedirs(output_dir,mode=0o755,exist_ok=True)
	# plt.savefig(f"{output_dir}/{file_prefix}.png",dpi=300)

	# return None

def run(args=None) -> None:

	plot_SE(SE_file=args.shannon_entropy_file,consensus_file=args.consensus_sequence_file,highlight_file=args.highlight_residues,output_dir=args.outdir)

	return None

if __name__ == "__main__":

	from argparse import ArgumentParser

	GetOptions = ArgumentParser()

	GetOptions.add_argument("-s","--shannon_entropy_file",required=True,type=str)
	GetOptions.add_argument("-o","--outdir",required=False,type=str,default="se_graphics")
	GetOptions.add_argument("-y","--highlight_residues",required=False,type=str)
	GetOptions.add_argument("-c","--consensus_sequence_file",required=False,type=str)

	run(GetOptions.parse_known_args()[0])