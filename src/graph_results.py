#!/usr/bin/env python3

import sys
import numpy as np
import pandas as pd
from math import log
from os import makedirs
from random import sample
from os.path import exists
from statistics import median, mean
import matplotlib.colors as mcolors
from matplotlib import pyplot as plt
import matplotlib.patches as patches
from matplotlib.ticker import MultipleLocator


def __read_configuration_file(config_file:str=None) -> dict:
	if not exists(config_file):
		print(f"\t\t[E]  Can't open configuration file located at {config_file}\n\n")
		return -1

	config = {}
	with open(config_file,'r') as IN:
		for line in IN:
			key,value=line.strip().split("=")
			config[key.strip()] = value.strip() == "True"

	return config


def __read_average_entropy_file(average_SE_file:str=None,msa_links:dict=None) -> np.array:
	
	if not exists(average_SE_file):
		print(f"\t\t[W]  Can't open average entropy file located at {average_SE_file}\n\n")
		return -1

	average_SEs = np.zeros(len(msa_links.keys()))
	with open(average_SE_file,'r') as IN:
		for line in IN:
			line = line.strip()
			if line == "":
				continue
			if line[0] == "#":
				continue
			res_num,average = line.split("\t")
			res_num = int(res_num)
			if res_num in msa_links.keys():
				average_SEs[msa_links[res_num]] = float(average)

	return average_SEs

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

def __read_res_highlight_subset_file(file:str=None) -> dict:

	if file is None or not exists(file):
		return None

	selection = {}

	source = None
	with open(file,'r') as IN:
		for line in IN:
			line = line.strip()

			if line == "":
				continue

			if line[0] == "#":
				source = line.replace("#","").strip()
				selection[source] = []
				continue

			values = line.split("-")

			if len(values) == 2:
				for x in range(int(values[0]),int(values[1])+1):
					if int(x) not in selection[source]:
						selection[source].append(int(x))
			else:
				if int(values[0]) not in selection[source]:
					selection[source].append(int(values[0]))

	return selection

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

			data[int(pos)] = {"RES":res,"SE":abs(float(fse)),"NUM_RES":None,"MSA_POS":msa_pos}
			msa_links[int(msa_pos)] = int(pos)

	return data,msa_links

def __read_zscale_data(zscale_file:str=None,data:dict=None,msa_links:dict=None) -> dict:

	if zscale_file is None or not exists(zscale_file):
		return None

	with open(zscale_file,'r') as IN:
		for line in IN:
			line = line.strip()

			if line == "":
				continue

			if line[0] == "#":
				continue

			msa_pos,z1,z2,z3,z4,z5 = map(float,line.split())
			msa_pos = int(msa_pos)

			if msa_pos in msa_links.keys():
				data[msa_links[msa_pos]]["ZSCALES"] = {"Z1":z1,"Z2":z2,"Z3":z3,"Z4":z4,"Z5":z5}
	
	return data

def read_data(SE_file:str=None,
			  consensus_file:str=None,
			  highlight_file:str=None,
			  average_SE_file:str=None,
			  subset_file:str=None,
			  zscale_file:str=None):
		
	# Read in Shannon Entropy
	data,msa_links = __read_SE_data(SE_file=SE_file)

	# Read in Average Shannon Entropy
	ase_values = __read_average_entropy_file(average_SE_file=average_SE_file,msa_links=msa_links)

	# Read in Consensus Sequence
	data = __read_consensus_file(consensus_file=consensus_file,data=data,msa_links=msa_links)

	# Read in residues to highlight in the graphic
	highlights = __read_res_highlight_subset_file(file=highlight_file)

	# Read in Z-Scale data
	zscales = __read_zscale_data(zscale_file=zscale_file,data=data,msa_links=msa_links)
	
	# Read subset file and filter data to only include subset residues
	subset = __read_res_highlight_subset_file(file=subset_file)
	if subset_file is not None:
		subset_pos = set(v for v in subset.values() for v in v)
		subset_indices = np.array(sorted(subset_pos), dtype=int)
		data = {key: value for key, value in data.items() if key in subset_pos}
		ase_values = ase_values[subset_indices]

	return data, ase_values, highlights
	

def plot(data:dict=None,
		 highlights:dict=None,
		 ase_values:np.array=None,
		 config_file=None,
		 fname_out=None) -> None:

	aas = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
	scales = ['Z1','Z2','Z3']

	# Read in configuration file
	config = __read_configuration_file(config_file=config_file)
	
	# Reference residue sequences
	res_nums = sorted(data.keys())

	# Number of sequences present in MSA
	num_of_seq = sum(x for x in data[res_nums[0]]["RES_COUNTS"].values())

	# Number of residues in the reference sequence
	len_of_seq = len(res_nums)

	num_ress = np.zeros(len_of_seq)

	# Font size of tick labels
	font_size = 6
	# Width of the graphic as to not have labels overlapping each other
	width = font_size*(1/72)*(len_of_seq)*1.5
	
	####################################################################################################################
	## Setting up graphics
	####################################################################################################################

	fig = plt.figure()

	# Determine the number of rows and columns for the plot from the configuration file
	plot_columns = 1
	
	# List of plots that take 1 row
	one_row_plots = ["se_graphing", "specificity_graphing", "number_of_residues", "residue_distribution", "zscale_distribution"]
	plot_rows = sum(1 for plot in one_row_plots if config[plot])

	# Add 2 rows if residue_presence_matrix is enabled
	if config["residue_presence_matrix"]:
		plot_rows += 2
	
	current_row = 0
	final_plot = None
	first_plot = None

	## Specifying plot placements
	# Shannon Entropy 1st
	if config["se_graphing"]:
		se_graph = plt.subplot2grid((plot_rows,plot_columns),(current_row,0))
		if current_row == 0:
			first_plot = se_graph
		current_row += 1
		final_plot = se_graph
	# Specificty/Conservation Scoring
	if config["specificity_graphing"]:
		scs_graph = plt.subplot2grid((plot_rows,plot_columns),(current_row,0))
		if current_row == 0:
			first_plot = scs_graph
		current_row += 1
		final_plot = scs_graph
	# Residue presense matrix 2nd, spans 2 rows
	if config["residue_presence_matrix"]:
		res_mat_graph = plt.subplot2grid((plot_rows,plot_columns),(current_row,0),rowspan=2)
		if current_row == 0:
			first_plot = res_mat_graph
		current_row += 2
		final_plot = res_mat_graph
	# Number of different residues present at site
	if config["number_of_residues"]:
		res_num_graph = plt.subplot2grid((plot_rows,plot_columns),(current_row,0))
		if current_row == 0:
			first_plot = res_num_graph
		current_row += 1
		final_plot = res_num_graph
	# Number of sequences belonging to that residue
	if config["residue_distribution"]:
		res_dist_graph = plt.subplot2grid((plot_rows,plot_columns),(current_row,0))
		if current_row == 0:
			first_plot = res_dist_graph
		current_row += 1
		final_plot = res_dist_graph
	# Z-scale distribution
	if config["zscale_distribution"]:
		zscale_graph = plt.subplot2grid((plot_rows,plot_columns),(current_row,0))
		if current_row == 0:
			first_plot = zscale_graph
		current_row += 1
		final_plot = zscale_graph

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
	se_values = np.array([data[x]["SE"] for x in res_nums])

	# Residue and residue number labels
	# labels = [f"{data[x]["RES"]} [{x+1}]" if x%2 != 0 else f"{data[x]["RES"]}" for x in res_nums]
	labels = [f"{data[x]['RES']} [{x+1}]" if x%2 != 0 else f"{data[x]['RES']}" for x in res_nums]

	# Highlight residues with near equal representation in MSA
	for i,val_x in enumerate(res_nums):
		num_res = data[val_x]["NUM_RES"]
		num_res -= 1 if '-' in data[val_x]["RES_COUNTS"].keys() else 0
		num_res -= 1 if 'X' in data[val_x]["RES_COUNTS"].keys() else 0
		num_ress[i] = num_res
		if num_res < 10:
			continue
		if config["se_graphing"]:
			se_graph.plot(val_x,data[val_x]["SE"],marker="*",markersize=10*(1-data[val_x]["DEV"]),color='b')

	if config["se_graphing"]:
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

		## Plot the data
		se_graph.plot(res_nums,se_values,linewidth=0.5,color='k')

		## Prettify the Y-axis
		se_graph.yaxis.set_ticks([i*0.2 for i in range(int(1/0.2)+1)])
		se_graph.yaxis.set_ticklabels([0.0,0.2,0.4,0.6,0.8,1.0])

		se_graph.tick_params(axis='x',which='both',bottom=False,labelbottom=False)
		se_graph.xaxis.grid(color='w',linestyle='-',linewidth=0.75,which='minor')
		se_graph.spines['bottom'].set_visible(False)

		# ase_graph = se_graph.twinx()
		se_graph.plot(res_nums,[ase_values[x] for x in res_nums],linewidth=0.5,color='b')

	# ####################################################################################################################
	# ## Specificity/Conservation Scoring
	# ####################################################################################################################

	rotate_deg = np.deg2rad(45)
	sin_val = np.sin(rotate_deg)
	cos_val = np.cos(rotate_deg)
	off = 0.5
	ase_values = ase_values
	h = np.sqrt(np.square(ase_values)+np.square(se_values))
	cons = np.zeros(h.shape)
	np.seterr(divide='ignore',invalid='ignore')
	cons = -h*np.cos(np.arccos(ase_values/h)-rotate_deg)
	np.nan_to_num(cons,copy=False)
	cons += off
	spec_dist = np.sqrt(np.square(ase_values)+np.square(np.subtract(se_values,1)))

	for index,x in enumerate(res_nums):
		# scs_graph.plot(x,cons[index],markersize=5*(1-spec_dist[index]),marker='s',color='orange')
		start_y = cons[index] if cons[index] < 0 else 0
		length_y = abs(cons[index])
		rect_width = 1-spec_dist[index]
		rect_start = x - (0.5*rect_width)
		if config["specificity_graphing"]:
			scs_graph.add_patch(plt.Rectangle((rect_start,start_y),rect_width,length_y,facecolor='b',edgecolor='b'))

	# ###########################################################
	# ## Hightlight user-provided residues
	# ###########################################################
	if highlights is not None:
		rand_colors = sample(range(len(color_keys)),len(highlights.keys()))
		for index,source in enumerate(highlights.keys()):
			for val_x in highlights[source]:
				if config["specificity_graphing"]:
					scs_graph.add_patch(plt.Rectangle((val_x-0.5,-0.6),1,1.2,facecolor=colors[f'xkcd:{color_keys[rand_colors[index]]}'],edgecolor="None",zorder=0,alpha=0.25))
	
	if config["specificity_graphing"]:
		# Set and rotate major tick labels
		if scs_graph == first_plot:
			labels = [None for val_x in res_nums]
		else:
			labels = [f"{val_x+1: >{len(str(len_of_seq))+1}}" if val_x%2 != 0 else "" for val_x in res_nums]
		scs_graph.xaxis.set_ticks([i for i in res_nums])
		scs_graph.tick_params(axis='x',which='both',bottom=False,labelbottom=False,top=False,labeltop=False)
		scs_graph.set_xticklabels(labels=labels,rotation='vertical',font='monospace',fontsize=font_size,va='center',ha='center')
		scs_graph.tick_params(axis='x',which='major',top=False,bottom=False,labelbottom=False,labeltop=True)
		scs_graph.xaxis.grid(color='w',linestyle='-',linewidth=0.75,which='minor')
		scs_graph.spines['bottom'].set_visible(False)
		scs_graph.spines['top'].set_visible(False)
		scs_graph.xaxis.set_minor_locator(MultipleLocator(1))
		scs_graph.grid(color='gainsboro',linestyle='-',linewidth=0.75,which='minor')
		scs_graph.set_xticks(res_nums)
		scs_graph.set_xlim([(res_nums[0]-0.5),res_nums[-1]+0.5])
		scs_graph.set_ylim([-0.65,0.65])

	# ####################################################################################################################
	# ## Residue Presence Matrix
	# ####################################################################################################################

	## Get a positive/negative matrix of what residues are present in the MSA at each residue
	present_res = np.zeros((len(aas),len(res_nums)))
	for row_index,aa in enumerate(aas):
		for col_index,val_x in enumerate(res_nums):
			if aa in data[val_x]["RES_COUNTS"].keys():
				present_res[row_index][col_index] = 1 if data[val_x]["RES_COUNTS"][aa] > 0 else 0

	if config["residue_presence_matrix"]:
		res_mat_graph.matshow(present_res,aspect='auto',cmap='tab20_r')

		## Setup Y-Axis
		res_mat_graph.set_yticks([i for i in range(len(aas))],labels=aas,font='monospace')
		# res_mat_graph.yaxis.set_minor_locator(MultipleLocator(1,offset=-0.5))
		res_mat_graph.yaxis.set_minor_locator(MultipleLocator(1))
		res_mat_graph.yaxis.grid(color='w',linestyle='-',linewidth=0.75,which='minor')
		res_mat_graph.tick_params(axis='y',which='minor',left=False)

		## Setup X-Axis
		# Set major tick for each residue [0 to pos-1]
		res_mat_graph.xaxis.set_ticks([i-1 for i in range(len(res_nums)+1)])
		# Set and rotate major tick labels
		if res_mat_graph == first_plot:
			labels = [None for val_x in res_nums]
		else:
			labels = [f"{val_x: >{len(str(len_of_seq))+1}}" if val_x%2 == 0 else "" for val_x in range(len(res_nums)+1)]
		res_mat_graph.set_xticklabels(labels=labels,rotation='vertical',font='monospace',fontsize=font_size,va='center',ha='center')
		# Show major tick labels, but not the ticks themselves, and only on the top
		res_mat_graph.tick_params(axis='x',which='major',top=False,bottom=False,labeltop=True,labelbottom=False)

		# Set minor ticks for each residue, but offset by -0.5 to create fake gridlines
		# res_mat_graph.xaxis.set_minor_locator(MultipleLocator(1,offset=-0.5))
		res_mat_graph.xaxis.set_minor_locator(MultipleLocator(1))
		# Don't show any information regarding the minor ticks
		res_mat_graph.tick_params(axis='x',which='minor',top=False,bottom=False,labeltop=False,labelbottom=False)
		
		# Hide the top and bottom axes
		res_mat_graph.spines['top'].set_visible(False)
		res_mat_graph.spines['bottom'].set_visible(False)

		# Turn on X-axis gridlines using the minor ticks
		res_mat_graph.xaxis.grid(color='w',linestyle='-',linewidth=0.75,which='minor')
		res_mat_graph.set_xlim([(res_nums[0]-0.5),res_nums[-1]+0.5])

	# ####################################################################################################################
	# ## Number of Residue Bar Graph
	# ####################################################################################################################

	if config["number_of_residues"]:
		res_num_graph.bar(res_nums,num_ress,1)

		## Setup X-Axis
		# Set major tick for each residue [0 to pos-1]
		res_num_graph.xaxis.set_ticks([i for i in res_nums])
		# Set and rotate major tick labels
		if res_num_graph == first_plot:
			labels = [None for val_x in res_nums]
		else:
			labels = [f"{val_x+1: >{len(str(len_of_seq))+1}}" if val_x%2 != 0 else "" for val_x in res_nums]
		res_num_graph.set_xticklabels(labels=labels,rotation='vertical',font='monospace',fontsize=font_size,va='center',ha='center')
		res_num_graph.tick_params(axis='x',which='major',top=False,bottom=False,labelbottom=False,labeltop=True)

		# res_num_graph.xaxis.set_minor_locator(MultipleLocator(1,offset=-0.5))
		res_num_graph.xaxis.set_minor_locator(MultipleLocator(1))
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
	
	# # ###################################################################################################################
	# # ## Residue Distribution Bar Graph
	# # ###################################################################################################################

	# TODO Add colorscheme here for residue classification (pos, neg, etc.)
	# Stacking bars require a y-value to place the bottom of the bar
	
	bottoms = np.zeros(len_of_seq)

	if config["residue_distribution"]:
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
		# res_dist_graph.xaxis.set_minor_locator(MultipleLocator(1,offset=-0.5))
		res_dist_graph.xaxis.set_minor_locator(MultipleLocator(1))
		res_dist_graph.grid(color='gainsboro',linestyle='-',linewidth=0.75,which='minor')

		## Setup Y-Axis
		res_dist_graph.yaxis.set_ticks([i*0.2 for i in range(int(1/0.2)+1)])
		res_dist_graph.yaxis.set_ticklabels([0.0,0.2,0.4,0.6,0.8,1.0])

		res_dist_graph.spines['top'].set_visible(False)

	# # ###################################################################################################################
	# # ## ZSCALES
	# # ###################################################################################################################

	if config["zscale_distribution"]:

		# Store Z-sclales in matrix
		m = np.zeros((len(scales), len(res_nums) +1))  # For clarity, show only first 3 Z-scales

		for i,scale in enumerate(scales):
			for index,val_x in enumerate(res_nums):
				m[i][index] = data[val_x]["ZSCALES"][scale]

		zscale_graph.matshow(m,aspect='auto',cmap='Reds')

		## Setup Y-Axis
		zscale_graph.set_yticks([i for i in range(len(scales))],labels=scales[::-1],font='monospace')

		## Setup X-Axis
		labels = [f"{val_x+1: >{len(str(len(res_nums)))+1}}" if val_x%2 != 0 else "" for val_x in res_nums]
		zscale_graph.set_xlim([(res_nums[0]-0.5),res_nums[-1]+0.5])
		# Set major tick for each residue [0 to pos-1]
		zscale_graph.xaxis.set_ticks([i for i in res_nums])
		
		# Set and rotate major tick labels
		zscale_graph.set_xticklabels(labels=labels,rotation='vertical',font='monospace',fontsize=font_size,va='center',ha='center')
		zscale_graph.tick_params(axis='x',which='major',top=False,bottom=False,labeltop=False,labelbottom=False, pad=10)
		zscale_graph.tick_params(axis='x',which='minor',top=False,bottom=False,labeltop=False,labelbottom=False, pad=10)

		zscale_graph.spines['bottom'].set_visible(False)
		zscale_graph.spines['top'].set_visible(False)
		zscale_graph.spines['left'].set_visible(False)
		zscale_graph.spines['right'].set_visible(False)

	## Setup top and bottom X-Axis
	# first_plot.set_xticks([])
	axt = first_plot.secondary_xaxis('top')
	axt.set_xticks(res_nums)
	axt.tick_params(which='both',top=False)
	labels = [f"{data[x]['RES']} [{x+1}]" if x%2 != 0 else f"{data[x]['RES']}" for x in res_nums]
	axt.set_xticklabels(labels,font='monospace',fontsize=font_size,ha='center',va='bottom',rotation='vertical')

	rdg_bt = final_plot.secondary_xaxis("bottom")
	rdg_bt.set_xticks([i for i in res_nums])
	rdg_bt.tick_params(which='both',bottom=False)
	# labels = [f"[{x+1}] {data[x]["CON_RES"]}" if x%2 != 0 else data[x]["CON_RES"] for x in res_nums]
	labels = [f"[{x+1}] {data[x]['CON_RES']}" if x%2 != 0 else data[x]['CON_RES'] for x in res_nums]
	rdg_bt.set_xticklabels(labels=labels,rotation='vertical',font='monospace',fontsize=font_size)

	# # ###################################################################################################################
	# # ## SAVING FILES
	# # ###################################################################################################################

	fig.set_size_inches(width,10)
	plt.tight_layout()
	plt.savefig(fname_out,dpi=500)
	
	return None

def run(args=None) -> None:

	ref = args.shannon_entropy_file.split("/")[-1].split(".")[0]
	makedirs(args.outdir,mode=0o755,exist_ok=True)

	data, ase_values, highlights = read_data(SE_file=args.shannon_entropy_file,
										     consensus_file=args.consensus_sequence_file,
											 highlight_file=args.highlight_residues,
											 average_SE_file=args.average_entropy_file,
											 subset_file=args.subset_file,
											 zscale_file=args.zscale_file)
	plot(data=data,
	     highlights=highlights,
		 ase_values=ase_values,
		 config_file=args.configuration,
		 fname_out=f"{args.outdir}/{ref}.png")
	
	return None

if __name__ == "__main__":

	from argparse import ArgumentParser

	GetOptions = ArgumentParser()

	GetOptions.add_argument("-s","--shannon_entropy_file",required=True,type=str)
	GetOptions.add_argument("-o","--outdir",required=False,type=str,default="se_graphics")
	GetOptions.add_argument("-y","--highlight_residues",required=False,type=str)
	GetOptions.add_argument("-c","--consensus_sequence_file",required=False,type=str)
	GetOptions.add_argument("-a","--average_entropy_file",required=False)
	GetOptions.add_argument("-z","--zscale_file",required=False,type=str)
	GetOptions.add_argument("-f","--subset_file",required=False,type=str)
	GetOptions.add_argument("-l","--configuration",required=False,type=str)

	run(args=GetOptions.parse_known_args()[0])