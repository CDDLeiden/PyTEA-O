import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

from PyTEAO.analysis.twoentropyanalysis import TwoEntropyAnalysis
from PyTEAO.visualization.subplotbase import SubplotBase
from PyTEAO.utils.visualization import register

@register("ShannonEntropyLine")
class ShannonEntropyPlot(SubplotBase):

	ROWSPAN = 1
	COLSPAN = 1

	def plot(self,axes:plt.axes):

		global_entropy_values = self.tea.global_entropy[self.reference_indicies]
		average_entropy_values = self.tea.average_entropy[self.reference_indicies]
		
		axes.plot(self.residue_numbers,[0.25 for x in self.residue_numbers],color='r',alpha=0.25)
		axes.plot(self.residue_numbers,[0.5 for x in self.residue_numbers],color='y',alpha=0.25)
		axes.plot(self.residue_numbers,[0.75 for x in self.residue_numbers],color='g',alpha=0.25)

		axes.plot(self.residue_numbers,global_entropy_values,linewidth=0.5,color=self.COLORS['maroon'],linestyle="-",marker='.',label='Global Entropy')
		axes.plot(self.residue_numbers,average_entropy_values,linewidth=0.5,color=self.COLORS['cyan'],linestyle="-",marker='.',label='Average Entropy')
		
		axes.set_ylim([-0.05,1.05])

		axes.yaxis.set_ticks([i*0.2 for i in range(int(1/0.2)+1)])
		axes.yaxis.set_ticklabels([0.0,0.2,0.4,0.6,0.8,1.0])

		axes.vlines([x+0.5 for x in range(len(self.residue_numbers)-1)],colors="#f9f9f9ff",ymin=0,ymax=1)

		axes.legend(loc='center right',bbox_to_anchor=(-0.01,0.5))

		super()._SubplotBase__setup_axes(axes)