import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from PyTEAO.analysis.twoentropyanalysis import TwoEntropyAnalysis
from PyTEAO.visualization.subplotbase import SubplotBase
from PyTEAO.utils.visualization import register

@register("ResiduePresenceMatrix")
class ResiduePresencePlot(SubplotBase):

	ROWSPAN = 2
	COLSPAN = 1

	def plot(self,axes:plt.axes):

		presence_matrix:pd.DataFrame = (self.msa.residue_counts.loc[self.reference_indicies,self.AAs]>0).astype("Sparse[int]")

		cmap_mask = np.array([i+1 for i in range(len(presence_matrix.columns))])

		gradient_presence_matrix = presence_matrix*cmap_mask[np.newaxis,:]

		axes.matshow(gradient_presence_matrix.T,aspect='auto',cmap=self.cmap)

		axes.set_yticks([i for i in range(len(presence_matrix.columns))],labels=self.AAs,font='monospace')

		axes.vlines([x+0.5 for x in range(len(self.reference_indicies)-1)],colors=self.GRID_GREY,ymin=-0.5,ymax=len(self.AAs)-0.5)
		axes.hlines([x+0.5 for x in range(len(self.AAs)-1)],colors=self.GRID_GREY,xmin=-0.5,xmax=len(self.reference_indicies)-0.5)

		axes.set_zorder(1_000_000)

		super()._SubplotBase__setup_axes(axes,offset=False)