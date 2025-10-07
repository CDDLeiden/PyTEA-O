import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from PyTEAO.analysis.twoentropyanalysis import TwoEntropyAnalysis
from PyTEAO.visualization.subplotbase import SubplotBase
from PyTEAO.utils.visualization import register

@register("ResidueDistributionBar")
class ResidueDistributionPlot(SubplotBase):

	ROWSPAN = 1
	COLSPAN = 1

	def plot(self,axes:plt.axes):

		bottoms = np.zeros(len(self.residue_numbers))

		color_keys = list(self.COLORS.keys())

		for i,aa in enumerate(self.AAs[::-1]):
			
			residue_fractions = self.msa.residue_counts.loc[self.reference_indicies][aa]/self.msa.num_sequences
			
			axes.bar(self.residue_numbers,residue_fractions,1,bottom=bottoms,color=self.COLORS[color_keys[-1-i]])
			bottoms += residue_fractions

		axes.yaxis.set_ticks([i*0.2 for i in range(int(1/0.2)+1)])
		axes.yaxis.set_ticklabels([0.0,0.2,0.4,0.6,0.8,1.0])

		axes.set_ylabel("Residue Distribution")

		axes.set_zorder(1_000_000)

		super()._SubplotBase__setup_axes(axes)