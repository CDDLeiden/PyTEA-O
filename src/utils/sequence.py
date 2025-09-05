

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

	def __init__(self,codon_definition_file:str=None):

		return None