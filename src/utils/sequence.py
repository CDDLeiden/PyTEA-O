

class SequenceUtilities:

	AAs = ['A','B','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','X','Z','-']

	aa_to_codon = {
		'A': ['GCT', 'GCC', 'GCA', 'GCG'],       # Alanine
		'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],  # Arginine
		'N': ['AAT', 'AAC'],                    # Asparagine
		'D': ['GAT', 'GAC'],                    # Aspartic acid
		'C': ['TGT', 'TGC'],                    # Cysteine
		'Q': ['CAA', 'CAG'],                    # Glutamine
		'E': ['GAA', 'GAG'],                    # Glutamic acid
		'G': ['GGT', 'GGC', 'GGA', 'GGG'],      # Glycine
		'H': ['CAT', 'CAC'],                    # Histidine
		'I': ['ATT', 'ATC', 'ATA'],             # Isoleucine
		'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],  # Leucine
		'K': ['AAA', 'AAG'],                    # Lysine
		'M': ['ATG'],                           # Methionine (Start)
		'F': ['TTT', 'TTC'],                    # Phenylalanine
		'P': ['CCT', 'CCC', 'CCA', 'CCG'],      # Proline
		'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],  # Serine
		'T': ['ACT', 'ACC', 'ACA', 'ACG'],      # Threonine
		'W': ['TGG'],                           # Tryptophan
		'Y': ['TAT', 'TAC'],                    # Tyrosine
		'V': ['GTT', 'GTC', 'GTA', 'GTG'],      # Valine
		'*': ['TAA', 'TAG', 'TGA'],             # Stop codons
	}

	codon_to_aa = {
		'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
		'TGT': 'C', 'TGC': 'C',
		'GAT': 'D', 'GAC': 'D',
		'GAA': 'E', 'GAG': 'E',
		'TTT': 'F', 'TTC': 'F',
		'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
		'CAT': 'H', 'CAC': 'H',
		'ATT': 'I', 'ATC': 'I', 'ATA': 'I',
		'TTA': 'L', 'TTG': 'L','CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
		'AAA': 'K', 'AAG': 'K',
		'ATG': 'M',
		'AAT': 'N', 'AAC': 'N',
		'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
		'CAA': 'Q', 'CAG': 'Q',
		'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R','AGA': 'R', 'AGG': 'R',
		'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
		'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S','AGT': 'S', 'AGC': 'S',
		'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
		'TAT': 'Y', 'TAC': 'Y',
		'TGG': 'W',
		'TAA': '*', 'TAG': '*', 'TGA': '*',
	}

	Sandberg_Zscales = {
		'A': {'Zscale_1': 0.24, 'Zscale_2': -2.32, 'Zscale_3': 0.6},
		'C': {'Zscale_1': 0.84, 'Zscale_2': -1.67, 'Zscale_3': 3.71},
		'D': {'Zscale_1': 3.98, 'Zscale_2': 0.93, 'Zscale_3': 1.93},
		'E': {'Zscale_1': 3.11, 'Zscale_2': 0.26, 'Zscale_3': -0.11},
		'F': {'Zscale_1': -4.22, 'Zscale_2': 1.94, 'Zscale_3': 1.06},
		'G': {'Zscale_1': 2.05, 'Zscale_2': -4.06, 'Zscale_3': 0.36},
		'H': {'Zscale_1': 2.47, 'Zscale_2': 1.95, 'Zscale_3': 0.26},
		'I': {'Zscale_1': -3.89, 'Zscale_2': -1.73, 'Zscale_3': -1.71},
		'K': {'Zscale_1': 2.29, 'Zscale_2': 0.89, 'Zscale_3': -2.49},
		'L': {'Zscale_1': -4.28, 'Zscale_2': -1.3, 'Zscale_3': -1.49},
		'M': {'Zscale_1': -2.85, 'Zscale_2': -0.22, 'Zscale_3': 0.47},
		'N': {'Zscale_1': 3.05, 'Zscale_2': 1.62, 'Zscale_3': 1.04},
		'P': {'Zscale_1': -1.66, 'Zscale_2': 0.27, 'Zscale_3': 1.84},
		'Q': {'Zscale_1': 1.75, 'Zscale_2': 0.5, 'Zscale_3': -1.44},
		'R': {'Zscale_1': 3.52, 'Zscale_2': 2.5, 'Zscale_3': -3.5},
		'S': {'Zscale_1': 2.39, 'Zscale_2': -1.07, 'Zscale_3': 1.15},
		'T': {'Zscale_1': 0.75, 'Zscale_2': -2.18, 'Zscale_3': -1.12},
		'V': {'Zscale_1': -2.59, 'Zscale_2': -2.64, 'Zscale_3': -1.54},
		'W': {'Zscale_1': -4.36, 'Zscale_2': 3.94, 'Zscale_3': 0.59},
		'Y': {'Zscale_1': -2.54, 'Zscale_2': 2.44, 'Zscale_3': 0.43}
	}

	def __init__(self,codon_definition_file:str=None):

		return None