import matplotlib.pyplot as plt

from src.analysis.twoentropyanalysis import TwoEntropyAnalysis
from src.visualization.subplotbase import SubplotBase
from src.utils.visualization import register

@register("SpecificityLine")
class SpecificityPlot(SubplotBase):

	ROWSPAN = 1
	COLSPAN = 1

	def plot(self,axes:plt.axes):
		
		conserved_scores = self.tea.conservation_scores[self.reference_indicies]
		specificity_scores = self.tea.specificity_scores[self.reference_indicies]

		axes.plot(self.residue_numbers,conserved_scores,color=self.COLORS['teal'],linewidth=0.5,marker=".",label='Residue Conservation')
		axes.plot(self.residue_numbers,specificity_scores,color=self.COLORS['purple'],linewidth=0.5,marker=".",label='Residue Specificity')
		axes.legend(loc='center right',bbox_to_anchor=(-0.01,0.5))
		axes.hlines(0,colors="#000000aa",xmin=self.residue_numbers[0]-0.5,xmax=self.residue_numbers[-1]+0.5)

		axes.set_ylim([-0.65,0.65])

		axes.vlines([x+0.5 for x in self.residue_numbers],colors="#f9f9f9ff",ymin=-.75,ymax=.75)

		super()._SubplotBase__setup_axes(axes)