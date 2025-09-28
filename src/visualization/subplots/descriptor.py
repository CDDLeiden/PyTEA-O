import matplotlib.pyplot as plt

from src.analysis.twoentropyanalysis import TwoEntropyAnalysis
from src.visualization.subplotbase import SubplotBase
from src.utils.visualization import register

@register("DescriptorHeatmap")
class DescriptorPlot(SubplotBase):

	ROWSPAN = 1
	COLSPAN = 1

	def plot(self,axes:plt.axes):


		deviation = self.msa.descriptor_deviation.iloc[self.reference_indicies].T
		axes.matshow(deviation,aspect='auto',cmap='Reds')

		ylabels = list(self.msa.descriptors_labels)

		## Setup Y-Axis
		axes.set_yticks([i for i in range(len(ylabels))],labels=ylabels[::-1],font='monospace')

		super()._SubplotBase__setup_axes(axes,offset=False)
		
		axes.set_zorder(1_000_000)