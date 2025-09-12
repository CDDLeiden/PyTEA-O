

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
		'A':[0.24,-2.32,0.6],
		'C':[0.84,-1.67,3.71],
		'D':[3.98,0.93,1.93],
		'E':[3.11,0.26,-0.11],
		'F':[-4.22,1.94,1.06],
		'G':[2.05,-4.06,0.36],
		'H':[2.47,1.95,0.26],
		'I':[-3.89,-1.73,-1.71],
		'K':[2.29,0.89,-2.49],
		'L':[-4.28,-1.3,-1.49],
		'M':[-2.85,-0.22,0.47],
		'N':[3.05,1.62,1.04],
		'P':[-1.66,0.27,1.84],
		'Q':[1.75,0.5,-1.44],
		'R':[3.52,2.5,-3.5],
		'S':[2.39,-1.07,1.15],
		'T':[0.75,-2.18,-1.12],
		'V':[-2.59,-2.64,-1.54],
		'W':[-4.36,3.94,0.59],
		'Y':[-2.54,2.44,0.43],
	}

	def __init__(self,codon_definition_file:str=None):

		return None