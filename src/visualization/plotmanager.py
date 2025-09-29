import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import typing
import pathlib
import random

from src.analysis.twoentropyanalysis import TwoEntropyAnalysis
from src.utils.visualization import import_all_subplots,SUBPLOT_REGISTRY
from src.visualization.subplotbase import SubplotBase
from src.utils.general import valid_directory

class PlotManager(SubplotBase):

	__SUBPLOTS = [
		"ShannonEntropyLine",
		"SpecificityLine",
		"ResiduePresenceMatrix",
		"UniqueResidueBar",
		"ResidueDistributionBar",
		"DescriptorHeatmap"
	]

	def __init__(self,tea:TwoEntropyAnalysis,subplots:str='111111',outdir:pathlib.Path|None="./plots",highlight_file:pathlib.Path|None=None):

		super().__init__(tea)

		import_all_subplots()

		self.tea = tea

		self.subplots = self.__decode_subplots_str(subplots)
		self.row,self.col = self.__get_gridspec_sizes()

		self.highlight_file = highlight_file

		self.outdir = valid_directory(outdir)

		self.width = self.FONT_SIZE*(1/72)*self.length_of_sequence*1.5

		self.fig = plt.figure(figsize=(self.width,self.row*2))

		self.gs = GridSpec(self.row,self.col,self.fig)

		current_row:int = 0
		self.axes:list[plt.axes] = []
		for plot in self.subplots:

			ax = self.fig.add_subplot(self.gs[current_row:current_row+plot.ROWSPAN,:])
			current_row += plot.ROWSPAN
			plot(self.tea).plot(ax)
			self.axes.append(ax)

		self.__add_top_labels(self.axes[0])
		self.__add_bottom_labels(self.axes[-1])

		if self.highlight_file is not None:
			self.__hightlight_residues(self.highlight_file)

		self.fig.tight_layout()


	def save_fig(self,file_type:str='svg',dpi:int=300):

		self.fig.savefig(self.outdir/f"figure.{file_type.replace(".","")}",dpi=dpi)

	def __get_highlight_residues(self,file:pathlib.Path) -> dict:

		selection = {}

		source:str
		with file.open('r') as IN:

			for line in IN:
				
				line = line.strip()

				if line == "":
					continue

				if line[0] == "#":
					source = line.replace("#","").strip()
					selection[source] = set()
					continue

				for value in line.split(","):
					split = value.split("-")
					try:
						split = [int(x) for x in split]
					except:
						raise TypeError(f"Specified residues must be specified by an integer, {split} is not valid.")
					span = list(range(split[0],split[-1]+1))
					selection[source].update(span)

		return selection

	def __hightlight_residues(self,file:pathlib.Path):

		highlight_set = self.__get_highlight_residues(file)

		### Pick random colors from the specified highlight colors
		rand_colors = random.sample(range(0,len(self.highlight_colors)-1),len(highlight_set.keys()))

		### Find the starting position and the width of the axes in the figure
		xmin,_,width,_ = self.axes[0].get_position(original=True).bounds

		### Calculate how wide each position in the MSA takes up in the figure
		step = width/(len(self.residue_numbers))


		for index,source in enumerate(highlight_set.keys()):

			for val_x in highlight_set[source]:

				## Skip out-of-bounds positions
				if val_x < 1 or val_x > self.length_of_sequence:
					continue

				## Create a residue stripe
				### Offset by +0.25 to center stripe
				### Offset by -1.0 to account for 0-indexing, residue 1 -> index 0
				### Stripe is only half the width of the residue
				residue_stripe = plt.Rectangle(
					(xmin+((val_x+0.25-1)*step),0),
					step/2,1,
					facecolor=self.highlight_colors[rand_colors[index]],
					edgecolor="None",
					zorder=0,
					alpha=0.25,
					transform=self.fig.transFigure
				)

				## Add residue stripe to the figure
				self.fig.patches.append(residue_stripe)

	def __add_top_labels(self,axes):

		ax = axes

		axt:plt.axes = axes.twiny()
		axt.set_xlim(ax.get_xlim())
		axt.set_xticks(ax.get_xticks())
		axt.tick_params(axis='x',which='major',top=False,bottom=False,labeltop=True,labelbottom=False)
		axt.tick_params(axis='x',which='minor',top=False,bottom=False,labeltop=False,labelbottom=False)
		
		axt.spines['bottom'].set_visible(False)
		axt.spines['top'].set_visible(True)

		top_labels = [f"{res} [{i+1}]" if i%2 != 0 else f"{res}" for i,res in enumerate(self.msa.get_sequence(self.reference_accesion))]

		axt.set_xticklabels(labels=top_labels,font='monospace',fontsize=self.FONT_SIZE,ha='center',va='bottom',rotation='vertical')
		axt.figure.canvas.draw()

	def __add_bottom_labels(self,axes):

		axb = axes
		consensus = self.msa.get_concensus_sequence()
		bottom_labels = [f"[{ref_index}] {consensus[ref_index]}" if i%2 != 0 else f"{consensus[ref_index]}" for i,ref_index in enumerate(self.reference_indicies)]
		axb.set_xticklabels(labels=bottom_labels,rotation='vertical',font='monospace',fontsize=self.FONT_SIZE,va='top',ha='center')
		axb.tick_params(axis='x',which='major',pad=5)
		axb.figure.canvas.draw()

	def __decode_subplots_str(self,subplots:str) -> list:

		__subplots = []

		for plot,setting in enumerate(subplots):

			if setting == "1":
				__subplots.append(SUBPLOT_REGISTRY.get(self.__SUBPLOTS[plot]))

		return __subplots

	def __get_gridspec_sizes(self) -> typing.Tuple[int,int]:

		row:int = 0
		col:int = 0

		for subplot in self.subplots:
			
			if subplot is None:
				continue
			
			row += getattr(subplot,"ROWSPAN",1)
			col += getattr(subplot,"COLSPAN",1)

		return row,col
