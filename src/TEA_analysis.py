#!/usr/bin/env python3

import os
import sys
import pickle
import numpy as np
import pandas as pd
import plotly.express as px
import matplotlib.pyplot as plt
from sklearn import preprocessing
from utils import extract_from_tsf, read_alignment, load_pickle

from matplotlib import pyplot as plt
from random import sample
import matplotlib.colors as mcolors
from matplotlib.ticker import MultipleLocator
import matplotlib.patches as patches
from os.path import exists
from os import makedirs
from statistics import median, mean
import numpy as np
import pandas as pd
from math import log
import sys
import matplotlib.gridspec as gridspec




class Analyse:
	"""
	Parse SE results and prepare for plotting
	"""

	def __init__(self):
		# TODO make some arguments necessary
		self.mode = None
		self.se_folder = None
		self.reference = None
		self.highlights_file = None

		self.SE_file = 'shannon_entropy.txt'
		self.average_SE_file = 'average_shannon_entropy.txt'
		self.summed_entropy_file = 'summed_subfamily_shannon_entropy.txt'
		self.consensus_file = 'consensus_logo.txt'
		self.zscale_file = 'zscales.txt'

		self.data = {}
		self.average_SEs: np.array = None
		self.summed_SEs: np.array = None
		self.highlights = None
		self.msa_links = {}

	def read_data(self):
		self.read_SE_data()
		self.read_consensus_file()
		self.read_zscale_file() # TODO make optional
	
		if self.mode == "TEA":
			self.summed_SEs = self.read_summed_entropy_file()

		elif self.mode == "TEAO":
			self.average_SEs = self.read_average_entropy_file()

		if self.highlights_file:
			self.highlights = self.read_res_highlight_file(self.highlights_file)
		
		# Save data to pickle
		with open(os.path.join(self.se_folder, f"{self.reference}_SE.pkl"), 'wb') as OUT:
			pickle.dump(self.data, OUT)


	def read_SE_data(self) -> list:

		SE_file = os.path.join(self.se_folder,  self.reference + '.' + self.SE_file)

		if not exists(SE_file):
			print(f"\n\n\t[E]  {SE_file} does not exist!\n")
			exit

		with open(SE_file,'r') as IN:
			for line in IN:
				line = line.strip()

				if line == "":
					continue

				if line[0] == "#":
					continue

				pos,msa_pos,res,se,fse,num_seqs,fog = line.split()

				self.data[int(pos)] = {"RES":res,"SE":abs(float(fse)),"NUM_RES":None,"MSA_POS":msa_pos}
				self.msa_links[int(msa_pos)] = int(pos)
		
		return self.data,self.msa_links


	def read_average_entropy_file(self) -> np.array:
		
		average_SE_file = os.path.join(self.se_folder, self.average_SE_file)

		if not exists(average_SE_file):
			print(f"\t\t[W]  Can't open average entropy file located at {average_SE_file}\n\n")
			return -1

		average_SEs = np.zeros(len(self.msa_links.keys()))
		with open(average_SE_file,'r') as IN:
			for line in IN:
				line = line.strip()
				if line == "":
					continue
				if line[0] == "#":
					continue
				res_num,average = line.split("\t")
				res_num = int(res_num)
				if res_num in self.msa_links.keys():
					average_SEs[self.msa_links[res_num]] = float(average)

		return average_SEs
	
	def read_summed_entropy_file(self) -> np.array:

		summed_entropy_file = os.path.join(self.se_folder, self.summed_entropy_file)

		if not exists(summed_entropy_file):
			print(f"\t\t[W]  Can't open summed entropy file located at {summed_entropy_file}\n\n")
			return -1
		
		summed_SEs = np.zeros(len(self.msa_links.keys()))
		with open(summed_entropy_file,'r') as IN:
			for line in IN:
				line = line.strip()
				if line == "":
					continue
				if line[0] == "#":
					continue
				res_num,summed = line.split("\t")
				res_num = int(res_num)
				if res_num in self.msa_links.keys():
					summed_SEs[self.msa_links[res_num]] = float(summed)

		return summed_SEs
	

	def read_consensus_file(self) -> list:

		consensus_file = os.path.join(self.se_folder, self.consensus_file)

		with open(consensus_file,'r') as IN:
			for line in IN:
				line = line.strip()

				if line == "":
					continue

				if line[0] == "#":
					continue

				msa_pos,res_nums = line.split()

				msa_pos = int(msa_pos)

				if msa_pos in self.msa_links.keys():
					non_gapped_data = [x for x in res_nums.split(";") if x.split(":")[0]]
					self.data[self.msa_links[msa_pos]]["CON_RES"] = non_gapped_data[0].split(":")[0]
					self.data[self.msa_links[msa_pos]]["RES_COUNTS"] = {x.split(":")[0]:float(x.split(":")[-1]) for x in non_gapped_data}
					self.data[self.msa_links[msa_pos]]["NUM_RES"] = len(non_gapped_data)
					temp_data = [int(x.split(":")[-1]) for x in non_gapped_data]
					self.data[self.msa_links[msa_pos]]["DEV"] = abs(mean(temp_data)-median(temp_data))/(mean(temp_data))
		
		return self.data

	def read_res_highlight_file(self, highlight_file:str=None) -> dict:

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
	
	def read_zscale_file(self) -> list:
	
		zscale_file = os.path.join(self.se_folder, self.zscale_file)

		with open(zscale_file,'r') as IN:
			for line in IN:
				line = line.strip()

				if line == "":
					continue

				if line[0] == "#":
					continue

				msa_pos,z1,z2,z3,z4,z5 = map(float,line.split())
				msa_pos = int(msa_pos)
	 
				if msa_pos in self.msa_links.keys():
					self.data[self.msa_links[msa_pos]]["ZSCALES"] = {"Z1":z1,"Z2":z2,"Z3":z3,"Z4":z4,"Z5":z5}
		
		return self.data


class Plot:
	"""
	Create some of the plots for the TEA analysis
	"""
	# TODO ADD TYPES HERE!! 
	def __init__(self, 
			  data, 
			  mode, 
			  subset=None,
			  summed_SEs=None,
			  average_SEs=None,
			  highlights=None):

		self.data = data
		self.mode = mode
		self.subset = subset
		self.summed_SEs = summed_SEs
		self.average_SEs = average_SEs
		self.highlights = highlights
		
		self.segment_threshold  = None
		self.segments = None
 
		self.res_nums = sorted(self.data.keys())  # Reference residue sequences

		# Filter data and only keep residues that are in the subset
		if self.subset:
			self.segment_threshold = 1.5
			self._filter_data()
		
		self.segments = self._segment_data()
		self.num_of_seq = sum(x for x in self.data[self.res_nums[0]]["RES_COUNTS"].values())  # Number of sequences present in MSA
		self.len_of_seq = len(self.res_nums)  # Number of residues in the reference sequence
		self.num_ress = np.zeros(self.len_of_seq)  # Number of different residues present at each site

		self.aas = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
		self.scales = ['Z1','Z2','Z3']

		# Font size of tick labels
		self.font_size = 6
		# Width of the graphic as to not have labels overlapping each other
		self.width = self.font_size*(1/72)*(self.len_of_seq)*1.5

		## Load CSS colors for iterative coloring
		self.colors = mcolors.XKCD_COLORS
		self.color_keys = [
			"aqua","blue","chartreuse","coral","crimson",
			"darkgreen","fuchsia","goldenrod","indigo","navy",
			"olive","orange","orangered","orchid","salmon",
			"teal","tomato","yellowgreen","grey","khaki"
		]

	def _filter_data(self):
		"""
		Only keep residues and residue-data that are in the provided subset.
		"""
		subset_values = set(value for values in self.subset.values() for value in values)
		self.res_nums = [resnum for resnum in self.res_nums if resnum in subset_values]
		# Loop over data and remove resnums that are not in subset
		if not self.average_SEs is None:
			tmp = [self.average_SEs[x] for x in self.res_nums]
			self.average_SEs = tmp
		if not self.summed_SEs is None:
			tmp = [self.summed_SEs[x] for x in self.res_nums]
			self.summed_SEs = tmp
		self.data = {k: v for k, v in self.data.items() if k in self.res_nums}
	

	def _segment_data(self):
		"""
		If subsets are provided, segments the data
		"""
		if self.segment_threshold is None:
			return [slice(0, len(self.res_nums))]  # Single segment covering all data

		# Identify gaps and segment the data
		gaps = np.diff(self.res_nums) > self.segment_threshold
		segment_indices = np.where(gaps)[0] + 1
		segments = []
		start_idx = 0

		for idx in segment_indices:
			segments.append(slice(start_idx, idx))
			start_idx = idx
		segments.append(slice(start_idx, len(self.res_nums)))

		return segments
	
	def _add_broken_axis_lines(self, axes):
		"""
		Adds diagonal lines to indicate breaks between subplots.
		
		Args:
			axes (list): List of axis objects.
		"""
		d = 0.015  # Size of diagonal lines in axis coordinates
		kwargs = dict(transform=axes[0].transAxes, color='k', clip_on=False)

		for i, ax in enumerate(axes):
			if i == 0:
				ax.plot((1 - d, 1 + d), (-d, +d), **kwargs)
				ax.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)
			elif i == len(axes) - 1:
				ax.plot((-d, +d), (1 - d, 1 + d), **kwargs)
				ax.plot((-d, +d), (-d, +d), **kwargs)
			else:
				kwargs.update(transform=ax.transAxes)
				ax.plot((-d, +d), (1 - d, 1 + d), **kwargs)
				ax.plot((-d, +d), (-d, +d), **kwargs)
				ax.plot((1 - d, 1 + d), (-d, +d), **kwargs)
				ax.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)


	def plot(self):

		####################################################################################################################
		## Setting up graphics
		####################################################################################################################

		master_fig = plt.figure(figsize=(30,20))

		n_segments = len(self.segments)

		plot_rows = 8
		plot_columns = n_segments

		# Correct height_ratios length to match plot_rows
		height_ratios = [1] * plot_rows
		width_ratios = [1] * plot_columns
		
		# Create the GridSpec layout for the subplots
		gs_master = gridspec.GridSpec(plot_rows, plot_columns, width_ratios=width_ratios, height_ratios=height_ratios, wspace=0.05, hspace=0.05)

		segment_index = 0
		for col in range(plot_columns):

			segment = self.segments[segment_index]

			# plot shannon entropy
			gs_pos = gs_master[0, col]
			ax = self.se_graphing(segment, figsize=(12, 12), se_graph=master_fig, gs_pos=gs_pos)

			# plot specificity
			gs_pos = gs_master[1, col]
			ax = self.specificity_graphing(segment, figsize=(12,12), scs_graph=master_fig, gs_pos=gs_pos)

			# self.highlight_residues(scs_graph)

			# plot residue presence matrix
			gs_pos = gs_master[2:4, col]
			ax = self.residue_presence_matrix(segment, figsize=(12,12), res_mat_graph=master_fig, gs_pos=gs_pos)

			# number of residues
			gs_pos = gs_master[5, col]
			ax = self.number_of_residues(segment, figsize=(12,12), res_num_graph=master_fig, gs_pos=gs_pos)

			# residue distribution
			gs_pos = gs_master[6, col]
			# ax = self.residue_distribution(segment, figsize=(12,12), res_dist_graph=master_fig, gs_pos=gs_pos)

			# zscale distribution
			gs_pos = gs_master[7, col]
			ax = self.zscale_distribution(segment, figsize=(12,12), zscale_graph=master_fig, gs_pos=gs_pos)

			segment_index += 1

		plt.tight_layout()


	def add_subplot(self, fig, gs_pos):
		ax = fig.add_subplot(gs_pos)
		ax.set_aspect('auto')
		return ax
	
	####################################################################################################################
	## Individual plotting functions
	####################################################################################################################

	def se_graphing(self, segment, figsize=(15, 5), se_graph=None, gs_pos=(0, 0)):

		# TODO Remove figsize and se_graph from function arguments

		# Shannon Entropies for the reference
		se_values = np.array([self.data[x]["SE"] for x in self.res_nums])

		ax = self.add_subplot(se_graph, gs_pos)

		if se_graph is None:
			se_graph = plt.figure(figsize=figsize)
			print("se_graphing, type se_graph",type(se_graph))

		# Plot
		x_val = self.res_nums[segment]
		indexes = [self.res_nums[segment].index(val) for val in x_val]
		y_val = [self.data[x]["SE"] for x in x_val]
		ax.plot(x_val, y_val, linewidth=0.5,color='k', zorder=10)
		y_val = [self.average_SEs[index] for index in indexes]
		ax.plot(x_val, y_val, linewidth=0.5,color='b')

		## Add "Highlighted Entropy Zones"
		# Green (0.0-0.25)
		# Yellow (0.25-0.75)
		# Red (0.75-1.0)
		ax.add_patch(plt.Rectangle((x_val[0]-1,-0.05),(x_val[-1]-x_val[0]+2),0.3,color='g',alpha=0.1))
		ax.add_patch(plt.Rectangle((x_val[0]-1,.25),(x_val[-1]-x_val[0]+2),0.5,color='y',alpha=0.1))
		ax.add_patch(plt.Rectangle((x_val[0]-1,.75),(x_val[-1]-x_val[0]+2),1,color='r',alpha=0.1))

		# Add separator lines to make residue visualization easier
		for x in x_val:
			ax.plot([x-0.5,x-0.5],[-0.05,1.05],color='w',linestyle='-',linewidth=0.75)

		# Highlight residues with near equal representation in MSA
		for i,val_x in enumerate(x_val):
			num_res = self.data[val_x]["NUM_RES"]
			num_res -= 1 if '-' in self.data[val_x]["RES_COUNTS"].keys() else 0
			num_res -= 1 if 'X' in self.data[val_x]["RES_COUNTS"].keys() else 0
			ax.plot(val_x,self.data[val_x]["SE"],marker="*",markersize=10*(1-self.data[val_x]["DEV"]),color='b')
			self.num_ress[i] = num_res
			if num_res < 10:
				continue

		## Prettify the Y-axis
		ax.yaxis.set_ticks([i*0.2 for i in range(int(1/0.2)+1)])
		ax.yaxis.set_ticklabels([0.0,0.2,0.4,0.6,0.8,1.0])

		ax.tick_params(axis='x',which='both',left=False,labelleft=False)
		ax.xaxis.grid(color='w',linestyle='-',linewidth=0.75,which='minor')
		ax.spines['bottom'].set_visible(False)

		ax.set_xlim([(x_val[0]-0.5),x_val[-1]+0.5])
		ax.set_ylim([-0.05,1.05])

	# 	# Residue and residue number labels
	# 	labels = [f"{self.data[x]['RES']} [{x+1}]" if x%2 != 0 else f"{self.data[x]['RES']}" for x in self.res_nums]

	# 	axt = se_graph.secondary_xaxis('top')
	# 	axt.set_xticks(self.res_nums)
	# 	axt.set_xticklabels(labels,font='monospace',fontsize=self.font_size,ha='center',va='bottom',rotation='vertical')
		
		# # Add vertical lines at breaks
		# for i in range(1, n_segments):  # Start from the second segment
		# 	# Get the rightmost x-coordinate of the previous segment
		# 	prev_segment = self.segments[i - 1]
		# 	prev_max_x = max(self.res_nums[prev_segment])  # Use the slice to get the max x
			
		# 	# Draw vertical line on the left of the current segment
		# 	axes[i].axvline(prev_max_x, color='black', linestyle='--', linewidth=2)

		# # Remove unnecessary ticks/labels
		# for i, ax in enumerate(axes):
		# 	if i > 0:
		# 		ax.set_ylabel('')  # Remove the y-label
		# 		ax.spines['left'].set_visible(False)
		# 		ax.tick_params(left=False, labelleft=False)

		return ax

	def specificity_graphing(self, segment, figsize, scs_graph, gs_pos):

		# Shannon Entropies for the reference
		se_values = np.array([self.data[x]["SE"] for x in self.res_nums])

		ax = self.add_subplot(scs_graph, gs_pos)

		x_val = self.res_nums[segment]
		indexes = [self.res_nums[segment].index(val) for val in x_val]

		rotate_deg = np.deg2rad(45)
		sin_val = np.sin(rotate_deg)
		cos_val = np.cos(rotate_deg)
		off = 0.5
		self.average_SEs = self.average_SEs
		h = np.sqrt(np.square(self.average_SEs)+np.square(se_values))
		cons = np.zeros(h.shape)
		np.seterr(divide='ignore',invalid='ignore')
		cons = -h*np.cos(np.arccos(self.average_SEs/h)-rotate_deg)
		np.nan_to_num(cons,copy=False)
		cons += off
		spec_dist = np.sqrt(np.square(self.average_SEs)+np.square(np.subtract(se_values,1)))

		for index,x in enumerate(x_val):
			ax.plot(x, cons[index],markersize=5*(1-spec_dist[index]),marker='s',color='orange')
			start_y = cons[index] if cons[index] < 0 else 0
			length_y = abs(cons[index])
			rect_width = 1-spec_dist[index]
			rect_start = x - (0.5*rect_width)
			ax.add_patch(plt.Rectangle((rect_start,start_y),rect_width,length_y,facecolor='b',edgecolor='b'))

		return ax

	def highlight_residues(self, scs_graph):

		if self.highlights:
			rand_colors = sample(range(len(self.color_keys)),len(self.highlights.keys()))
			for index,source in enumerate(self.highlights.keys()):
				for val_x in self.highlights[source]:
					scs_graph.add_patch(plt.Rectangle((val_x-0.5,-0.6),1,1.2,facecolor=self.colors[f'xkcd:{self.color_keys[rand_colors[index]]}'],edgecolor="None",zorder=0,alpha=0.25))

		# Set and rotate major tick labels
		labels = [f"{val_x+1: >{len(str(self.len_of_seq))+1}}" if val_x%2 != 0 else "" for val_x in self.res_nums]
		scs_graph.xaxis.set_ticks([i for i in self.res_nums])
		scs_graph.tick_params(axis='x',which='both',bottom=False,labelbottom=False,top=False,labeltop=False)
		scs_graph.set_xticklabels(labels=labels,rotation='vertical',font='monospace',fontsize=self.font_size,va='center',ha='center')
		scs_graph.tick_params(axis='x',which='major',top=False,bottom=False,labelbottom=False,labeltop=True)
		scs_graph.xaxis.grid(color='w',linestyle='-',linewidth=0.75,which='minor')
		scs_graph.spines['bottom'].set_visible(False)
		scs_graph.spines['top'].set_visible(False)
		scs_graph.xaxis.set_minor_locator(MultipleLocator(1))
		scs_graph.grid(color='gainsboro',linestyle='-',linewidth=0.75,which='minor')
		scs_graph.set_xticks(self.res_nums)
		scs_graph.set_xlim([(self.res_nums[0]-0.5),self.res_nums[-1]+0.5])
		scs_graph.set_ylim([-0.65,0.65])

	def residue_presence_matrix(self, segment, figsize, res_mat_graph, gs_pos):

		ax = self.add_subplot(res_mat_graph, gs_pos)

		x_val = self.res_nums[segment]
		indexes = [self.res_nums[segment].index(val) for val in x_val]

		aas = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']

		## Get a positive/negative matrix of what residues are present in the MSA at each residue
		present_res = np.zeros((len(aas),len(x_val)))
		for row_index,aa in enumerate(aas):
			for col_index,val_x in enumerate(x_val):
				if aa in self.data[val_x]["RES_COUNTS"].keys():
					present_res[row_index][col_index] = 1 if self.data[val_x]["RES_COUNTS"][aa] > 0 else 0

		# print(f"Present Residue Matrix: {present_res}")
		print(present_res.shape)
		print(present_res)


		ax.cla()
		img = ax.imshow(present_res,aspect='auto',cmap='tab20_r')
		img.set_data(present_res)

		## Setup Y-Axis
		ax.set_yticks([i for i in range(len(aas))],labels=aas,font='monospace')
		# res_mat_graph.yaxis.set_minor_locator(MultipleLocator(1,offset=-0.5))
		ax.yaxis.set_minor_locator(MultipleLocator(1))
		ax.yaxis.grid(color='w',linestyle='-',linewidth=0.75,which='minor')
		ax.tick_params(axis='y',which='minor',left=False)

		## Setup X-Axis
		# Set major tick for each residue [0 to pos-1]
		# ax.xaxis.set_ticks([i-1 for i in range(len(x_val)+1)])
		# Set and rotate major tick labels
		# labels = [f"{val_x: >{len(str(self.len_of_seq))+1}}" if val_x%2 == 0 else "" for val_x in range(len(x_val)+1)]
		# ax.set_xticklabels(labels=labels,rotation='vertical',font='monospace',fontsize=self.font_size,va='center',ha='center')
		# Show major tick labels, but not the ticks themselves, and only on the top
		ax.tick_params(axis='x',which='major',top=False,bottom=False,labeltop=True,labelbottom=False)

		# Set minor ticks for each residue, but offset by -0.5 to create fake gridlines
		# res_mat_graph.xaxis.set_minor_locator(MultipleLocator(1,offset=-0.5))
		ax.xaxis.set_minor_locator(MultipleLocator(1))
		# Don't show any information regarding the minor ticks
		ax.tick_params(axis='x',which='minor',top=False,bottom=False,labeltop=False,labelbottom=False)
		
		# Hide the top and bottom axes
		ax.spines['top'].set_visible(False)
		ax.spines['bottom'].set_visible(False)

		# Turn on X-axis gridlines using the minor ticks
		ax.xaxis.grid(color='w',linestyle='-',linewidth=0.75,which='minor')
		# ax.set_xlim([(x_val[0]-0.5),x_val[-1]+0.5])
		
		return ax

	def number_of_residues(self, segment, figsize, res_num_graph, gs_pos):
		"""
		Plot the number of different residues present at each position
		"""

		ax = self.add_subplot(res_num_graph, gs_pos)

		x_val = self.res_nums[segment]
		indexes = [self.res_nums[segment].index(val) for val in x_val]

		n_res = [self.num_ress[index] for index in indexes]
	
		ax.bar(x_val,n_res,1)

		## Setup X-Axis
		# Set major tick for each residue [0 to pos-1]
		ax.xaxis.set_ticks([i for i in x_val])
		# Set and rotate major tick labels
		labels = [f"{val_x+1: >{len(str(self.len_of_seq))+1}}" if val_x%2 != 0 else "" for val_x in x_val]
		ax.set_xticklabels(labels=labels,rotation='vertical',font='monospace',fontsize=self.font_size,va='center',ha='center')
		ax.tick_params(axis='x',which='major',top=False,bottom=False,labelbottom=False,labeltop=True)

		# res_num_graph.xaxis.set_minor_locator(MultipleLocator(1,offset=-0.5))
		ax.xaxis.set_minor_locator(MultipleLocator(1))
		ax.grid(color='gainsboro',linestyle='-',linewidth=0.75,which='minor')
		ax.tick_params(axis='x',which='minor',top=False,bottom=False)

		ax.set_xlim([(x_val[0]-0.5),x_val[-1]+0.5])

		## Setup Y-Axis
		ax.set_ylim([-0.05,21])
		ax.set_yticks([0,5,10,15,20])
		
		# Show horizontal guidelines
		ax.plot([x_val[0]-0.5,x_val[-1]+0.5],[5,5],color='w',linestyle='--',linewidth=1)
		ax.plot([x_val[0]-0.5,x_val[-1]+0.5],[10,10],color='aqua',linestyle='-',linewidth=1)
		ax.plot([x_val[0]-0.5,x_val[-1]+0.5],[15,15],color='w',linestyle='--',linewidth=1)
		ax.plot([x_val[0]-0.5,x_val[-1]+0.5],[20,20],color='aqua',linestyle='-',linewidth=1)

		# Hide the top and bottom axes
		ax.spines['top'].set_visible(False)
		ax.spines['bottom'].set_visible(False)

		return ax

	def residue_distribution(self, segment, figsize, res_dist_graph, gs_pos):

		# TODO ADD COLORSCHEME HERE

		# Stacking bars require a y-value to place the bottom of the bar

		aas = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']

		ax = self.add_subplot(res_dist_graph, gs_pos)

		x_val = self.res_nums[segment]
		indexes = [self.res_nums[segment].index(val) for val in x_val]
	
		bottoms = np.zeros(self.len_of_seq)
	
		labels = [f"{val_x+1: >{len(str(self.len_of_seq))+1}}" if val_x%2 != 0 else "" for val_x in x_val]

		for i in range(len(aas)):
			res_vals = np.zeros(self.len_of_seq)
			for index,val_x in enumerate(x_val):
				if i < len(self.data[val_x]["RES_COUNTS"].keys()):
					sorted_res_counts = sorted(self.data[val_x]["RES_COUNTS"].keys(),reverse=True,key=lambda x: self.data[val_x]["RES_COUNTS"][x])
					res_vals[index] = self.data[val_x]["RES_COUNTS"][sorted_res_counts[i]]/self.num_of_seq if sorted_res_counts[i] != "-" else 0
				else:
					res_vals[index] = 0

		# 	ax.bar(x_val,res_vals[],1,bottom=bottoms)
		# bottoms += res_vals

		## Setup top X-Axis
		ax.set_xlim([(x_val[0]-0.5),x_val[-1]+0.5])
		# Set major tick for each residue [0 to pos-1]
		ax.xaxis.set_ticks([i for i in x_val])
		# Set and rotate major tick labels
		ax.set_xticklabels(labels=labels,rotation='vertical',fontsize=self.font_size,va='center',ha='center')
		ax.tick_params(axis='x',which='major',top=False,bottom=False,labelbottom=False,labeltop=True)
		ax.tick_params(axis='x',which='minor',top=False,bottom=False,labeltop=False,labelbottom=False)
		# res_dist_graph.xaxis.set_minor_locator(MultipleLocator(1,offset=-0.5))
		ax.xaxis.set_minor_locator(MultipleLocator(1))
		ax.grid(color='gainsboro',linestyle='-',linewidth=0.75,which='minor')
		## Setup bottom X-Axis
		rdg_bt = ax.secondary_xaxis("bottom")
		rdg_bt.set_xticks([i for i in x_val])
		rdg_bt.tick_params(which='both',bottom=False)
		# labels = [f"[{x+1}] {data[x]["CON_RES"]}" if x%2 != 0 else data[x]["CON_RES"] for x in res_nums]
		labels = [f"[{x+1}] {self.data[x]['CON_RES']}" if x%2 != 0 else self.data[x]['CON_RES'] for x in x_val]
		rdg_bt.set_xticklabels(labels=labels,rotation='vertical',font='monospace',fontsize=self.font_size)

		## Setup Y-Axis
		ax.yaxis.set_ticks([i*0.2 for i in range(int(1/0.2)+1)])
		ax.yaxis.set_ticklabels([0.0,0.2,0.4,0.6,0.8,1.0])

		ax.spines['top'].set_visible(False)

		return ax

	def zscale_distribution(self, segment, figsize, zscale_graph, gs_pos):

		scales = ['Z1','Z2','Z3']

		ax = self.add_subplot(zscale_graph, gs_pos)

		x_val = self.res_nums[segment]
		indexes = [self.res_nums[segment].index(val) for val in x_val]

		# Store Z-sclales in matrix
		m = np.zeros((len(scales), len(x_val) +1))  # For clarity, show only first 3 Z-scales

		for i,scale in enumerate(scales):
			for index,val_x in enumerate(x_val):
				m[i][index] = self.data[val_x]["ZSCALES"][scale]

		print(f"Segment: {segment}, Data Shape: {m.shape}")

		# clear the axis
		ax.cla()
		img = ax.imshow(m, aspect='auto',cmap='Reds')
		img.set_data(m)

		## Setup Y-Axis
		ax.set_yticks([i for i in range(len(scales))],labels=scales[::-1],font='monospace')

		# ## Setup X-Axis
		ax.set_xticks(range(len(x_val)))
		ax.set_xticklabels(x_val, rotation=45, ha='right', font='monospace')

		# Set and rotate major tick labels
		# labels = [f"{val_x: >{len(str(self.len_of_seq))+1}}" if val_x%2 != 0 else "" for val_x in x_val]
		# ax.set_xticklabels(labels=labels,rotation='vertical',font='monospace',fontsize=self.font_size,va='center',ha='center')
		# Show major tick labels, but not the ticks themselves, and only on the top
		# ax.tick_params(axis='x',which='major',top=False,bottom=False,labeltop=False,labelbottom=True)

		# # Set minor ticks for each residue, but offset by -0.5 to create fake gridlines
		# res_mat_graph.xaxis.set_minor_locator(MultipleLocator(1,offset=-0.5))
		# ax.xaxis.set_minor_locator(MultipleLocator(1))
		# Don't show any information regarding the minor ticks
		# ax.tick_params(axis='x',which='minor',top=False,bottom=False,labeltop=False,labelbottom=True)
		
		# # Hide the top and bottom axes
		ax.spines['top'].set_visible(False)
		ax.spines['bottom'].set_visible(False)

		# # Turn on X-axis gridlines using the minor ticks
		ax.xaxis.grid(color='w',linestyle='-',linewidth=0.75,which='minor')

		return ax


	def global_vs_subfam(self, figsize=(5, 5)):
		f, ax = plt.subplots(figsize=figsize)
		if self.mode == "TEA":
			ax.scatter(self.summed_SEs, np.array([self.data[i]['SE'] for i in self.data]))
			ax.set_xlabel('Summed subfamily entropy')
		elif self.mode == "TEAO":
			ax.scatter(self.average_SEs, np.array([self.data[i]['SE'] for i in self.data]))
			ax.set_xlabel('Average subfamily entropy')
		ax.set_ylabel('Global entropy')
		ax.legend()
		return f, ax
	
	def global_vs_subfam_interactive(self, width=800, height=800):
		if self.mode == "TEAO":
			fig = px.scatter(x=self.average_SEs, 
							 y=np.array([self.data[i]['SE'] for i in self.data]),
							hover_name=np.array([self.data[i]['RES'] for i in self.data]),
							width=width,
							height=height)
			fig.update_layout(title='Average vs Global subfamily entropy', xaxis_title='Average subfamily entropy', yaxis_title='Global entropy')
		elif self.mode == "TEA":
			fig = px.scatter(x=self.summed_SEs, 
					y=np.array([self.data[i]['SE'] for i in self.data]),
					hover_name=np.array([self.data[i]['RES'] for i in self.data]),
					width=width,
					height=height)
			fig.update_layout(title='Summed vs Global subfamily entropy', xaxis_title='Summed subfamily entropy', yaxis_title='Global entropy')
		return fig
	
	# def plot_with_broken_axis(self, scs_graph=None, figsize=(15, 5), **kwargs):
	# 	"""
	# 	Creates a plot with broken axes if the data contains large gaps and adds vertical lines at breaks.

	# 	Args:
	# 		scs_graph (matplotlib.axes.Axes): Existing subplot, if any.
	# 		figsize (tuple): Figure size, e.g., (15, 5).
	# 		**kwargs: Additional plotting options.
	# 	"""
	# 	# Segment data
	# 	self.segments = self._segment_data()  # Assume this returns slices
	# 	n_segments = len(self.segments)

	# 	# Create figure using gridspec for better control
	# 	fig = plt.figure(figsize=figsize)
	# 	gs = gridspec.GridSpec(1, n_segments, width_ratios=[1] * n_segments, wspace=0.1)
	# 	axes = [fig.add_subplot(gs[i]) for i in range(n_segments)]

	# 	# Plot each segment
	# 	for ax, segment in zip(axes, self.segments):
	# 		self._plot_segment(ax, segment, **kwargs)

	# 	# Add vertical lines at breaks
	# 	for i in range(1, n_segments):  # Start from the second segment
	# 		# Get the rightmost x-coordinate of the previous segment
	# 		prev_segment = self.segments[i - 1]
	# 		prev_max_x = max(self.res_nums[prev_segment])  # Use the slice to get the max x

	# 		# Draw vertical line on the left of the current segment
	# 		axes[i].axvline(prev_max_x, color='black', linestyle='--', linewidth=2)

	# 	# Remove unnecessary ticks/labels
	# 	for i, ax in enumerate(axes):
	# 		if i > 0:
	# 			ax.set_ylabel('')  # Remove the y-label
	# 			ax.spines['left'].set_visible(False)
	# 			ax.tick_params(left=False, labelleft=False)

	# 	# Add broken axis lines
	# 	self._add_broken_axis_lines(axes)

	# 	# Final adjustments
	# 	plt.tight_layout(pad=0)