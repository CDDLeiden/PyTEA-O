import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from PyTEAO.analysis.twoentropyanalysis import TwoEntropyAnalysis
from PyTEAO.visualization.subplotbase import SubplotBase
from PyTEAO.utils.visualization import register

@register("UniqueResidueBar")
class UniqueResiduesPlot(SubplotBase):

	ROWSPAN = 1
	COLSPAN = 1

	def plot(self,axes:plt.axes):

		mask = ((self.msa.residue_counts>0).loc[self.reference_indicies][self.AAs]).astype("Sparse[int]")

		unique_counts = mask.sum(axis=1)
		
		axes.bar(self.residue_numbers,unique_counts.T,1,color=self.COLORS['navy'])

		# ## Setup Y-Axis
		axes.set_ylim([-0.05,20.5])
		axes.set_yticks([0,5,10,15,20])
		
		# # Show horizontal guidelines
		axes.hlines([5],colors=self.GRID_GREY,linestyle='--',xmin=0,xmax=len(self.residue_numbers))
		axes.hlines([10],colors=self.COLORS['cyan'],linestyle='-',xmin=0,xmax=len(self.residue_numbers),zorder=0)
		axes.hlines([15],colors=self.GRID_GREY,linestyle='--',xmin=0,xmax=len(self.residue_numbers))
		axes.hlines([20],colors=self.COLORS['cyan'],linestyle='-',xmin=0,xmax=len(self.residue_numbers),zorder=0)
		axes.set_ylabel("# of Unique Residues")

		# # Hide the top and bottom axes
		axes.spines['top'].set_visible(False)
		axes.spines['bottom'].set_visible(False)

		axes.vlines([x+0.5 for x in range(len(self.reference_indicies)-1)],colors=self.GRID_GREY,ymin=0,ymax=len(self.AAs))

		axes.set_zorder(1_000_000)
		super()._SubplotBase__setup_axes(axes)