import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

from src.utils.msa import MSA
from src.analysis.twoentropyanalysis import TwoEntropyAnalysis
from src.utils.sequence import SequenceUtilities


class SubplotBase:

	AAs = SequenceUtilities.Natural_AAs

	FONT_SIZE = 6

	COLORS = {
		"white":"#ffffffff",
		"maroon":"#800000",
		"brown":"#9A6324",
		"olive":"#808000",
		"teal":"#469990",
		"navy":"#000075",
		"red":"#e6194B",
		"orange":"#f58231",
		"yellow":"#ffe119",
		"lime":"#bfef45",
		"green":"#3cb44b",
		"cyan":"#42d4f4",
		"blue":"#4363d8",
		"purple":"#911eb4",
		"magenta":"#f032e6",
		"grey":"#a9a9a9",
		"pink":"#fabed4",
		"apricot":"#ffd8b1",
		"beige":"#fffac8",
		"mint":"#aaffc3",
		"lavender":"#dcbeff",
	}

	GRID_GREY = '#f9f9f9ff'

	def __init__(self,tea:TwoEntropyAnalysis):

		self.tea:TwoEntropyAnalysis = tea

		self.msa:MSA = tea.msa

		self.reference_accesion:str = self.msa.reference_accession

		self.reference_indicies = self.msa.get_reference_indices(self.reference_accesion).to_list()
		self.length_of_sequence = len(self.reference_indicies)
		self.residue_numbers = [x for x in range(self.length_of_sequence)]

		self.labels = [(x+1) if (x+1)%2 == 0 else "" for x in self.residue_numbers]

		self.__init_colors()


	def __init_colors(self):

		self.cmap = mcolors.ListedColormap([x for x in self.COLORS.values()])

		self.highlight_colors = [
			self.COLORS['maroon'],
			self.COLORS['brown'],
			self.COLORS['olive'],
			self.COLORS['teal'],
			self.COLORS['navy'],
			self.COLORS['red'],
			self.COLORS['orange'],
			self.COLORS['green'],
			self.COLORS['blue'],
			self.COLORS['purple'],
			self.COLORS['magenta'],
		]

	def __setup_axes(self,axes:plt.axes,offset:bool=True):

		## Setup X-Axis
		axes.set_xticks(self.residue_numbers)
		
		# # Set and rotate major tick labels
		axes.set_xticklabels(labels=self.labels,rotation='vertical',font='monospace',fontsize=self.FONT_SIZE,va='center',ha='center')
		axes.tick_params(axis='x',which='major',top=False,bottom=False,labeltop=False,labelbottom=True,length=0,pad=10)
		axes.tick_params(axis='x',which='minor',top=False,bottom=False,labeltop=False,labelbottom=False,length=0)
		axes.spines['bottom'].set_visible(False)
		axes.spines['top'].set_visible(False)

		if offset:

			axes.set_xlim([-0.5,self.residue_numbers[-1]+0.5])