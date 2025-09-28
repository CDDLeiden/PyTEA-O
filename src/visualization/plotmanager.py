import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import typing
import pathlib

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

	def __init__(self,tea:TwoEntropyAnalysis,subplots:str='111111',outdir:pathlib.Path|None="./plots"):

		super().__init__(tea)

		import_all_subplots()

		self.subplots = self.__decode_subplots_str(subplots)
		self.row,self.col = self.__get_gridspec_sizes()

		self.outdir = valid_directory(outdir)

		self.width = self.FONT_SIZE*(1/72)*self.length_of_sequence*1.5

		fig = plt.figure(figsize=(self.width,self.row*2))

		gs = GridSpec(self.row,self.col,fig)

		current_row:int = 0
		axes:list[plt.axes] = []
		for plot in self.subplots:

			ax = fig.add_subplot(gs[current_row:current_row+plot.ROWSPAN,:])
			current_row += plot.ROWSPAN
			plot(tea).plot(ax)
			axes.append(ax)

		self.__add_top_labels(axes[0])
		self.__add_bottom_labels(axes[-1])
		

		plt.tight_layout()
		plt.savefig(f"{outdir}/test.png",dpi=300)

	def __hightlight_residues(self):

		...
		# if highlights is not None:
		# rand_colors = sample(range(0,len(highlight_colors)-1),len(highlights.keys()))
		# xmin,ymax,width,height = plots[0].get_position(original=True).bounds
		# _,ymin,_,_ = plots[-1].get_position(original=True).bounds
		# step = width/(len(res_nums))
		# for index,source in enumerate(highlights.keys()):
		# 	for val_x in highlights[source]:
		# 		fig.patches.append(plt.Rectangle((xmin+((val_x+0.25)*step),0),step/2,1,facecolor=highlight_colors[rand_colors[index]],edgecolor="None",zorder=0,alpha=0.25,transform=fig.transFigure))

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
