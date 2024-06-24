import os
import sys
import math
import argparse
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio import Align
from Bio import AlignIO
from collections import Counter
from Bio.Align import substitution_matrices
from utils import extract_from_tsf, pairwise_seqidentity_matrix, blosum62_scores_matrix, save_pickle, rename, dict_to_json, save_commandline_args


class SE:
    """
    Calculate Shannon Entropy for a Multiple Sequence Alignment (MSA)

    """

    def __init__(self):

        # Command line arguments
        self.path_msa = None
        self.path_msa_tsv = None
        self.path_subfamilies = None  # Optional argument
        self.outpath = None  # It is important that the results are written in an empty directory. Make sure that the directory is empty and does no contain e.g. your alignment.
        self.headers = False
        self.reference_id = None
        self.calculate_pairwise_identity = False
        self.blosum = False

        # Default settings
        self.gap_symbol = "-"  # Gap symbol in MSA
        self.threshold_sequences_per_subfamily = 1  # Minimalsequences per subfamily
        self.shannon_entropy_file = "shannon_entropy.txt"
        self.shannon_entropy_file_filtered = "shannon_entropy_filtered.txt"
        self.shannon_entropy_file_header = ["alignment_pos", "entropy", "n_gaps", "n_nogap", "gap_fraction"]
        
        # For calculating similarity between subgroups, based on blosum62 matrix
        self.similarity_file = 'similarity_blosum62.txt'
        self.similarity_file_filtered = "similarity_blosum62_filtered.txt"
        self.similarity_file_header = ["alignment_pos", "blosum62_score", "n_gaps", "n_nogap", "gap_fraction"]

        self.n_res = None  # Ncols 
        self.n_seq = None  # Nrows

        self.msa_type = None

    def run(self):

        # Check if there already are entropy files in the output dir > otherwise ask to remove the present files or choose another directory
        if os.path.isfile(os.path.join(self.outpath, self.shannon_entropy_file_filtered)):
            sys.exit(f"Entropy file already exists in {self.outpath}, please choose another output directory or remove files.")

        # If the input is a .fasta file > convert to tab-separated file for calculations
        if SE.is_fasta(self.path_msa):
            self.msa_type = 'fasta'
            SE.save_tab_separated_msa_totxt(self.path_msa)
            self.path_msa_tsv = SE.get_tab_separated_fname(self.path_msa)

        # Assume that if input file is not a fasta file, it is a tab-separated file
        if not self.path_msa_tsv:
            self.msa_type = 'tab-separated'
            self.path_msa_tsv = self.path_msa
 
        # Make sure all sequences in the MSA have the same length before starting calculations
        result, length_seqs = self.check_size_msa()  
        if not result:
            sys.exit("MSA does not have the same length for all sequences")

        # Set matrix size based on size of input MSA
        self.set_size_msa(length_seqs)

        # Make dataframe with a separate column for every residue position in the MSA
        M, ids = self.format_data()  
        df = self.make_df(M, ids)

        # Create output file and write header to it
        SE_outputfile = os.path.join(self.outpath, self.shannon_entropy_file)
        with open(SE_outputfile, 'w+') as f:
            f.write('\t'.join(self.shannon_entropy_file_header) + '\n')

        if self.blosum:
            blosum_outputfile = os.path.join(self.outpath, self.similarity_file)
            with open(blosum_outputfile, 'w+') as f:
                f.write('\t'.join(self.similarity_file_header) + '\n')

        # Calculate Shannon Entropy per residue position (global entropy)
        for i in range(1, self.n_res + 1):
            col = f'Pos_{i}'  # Starting from 1
            arr = df[col].tolist()
            n_gaps = arr.count(self.gap_symbol)  # Count and remove gaps
            gap_fraction = round(n_gaps / len(arr), 2)
            # filtered = list(filter(lambda x: x != self.gap_symbol, arr))
            filtered = self.replace_gaps(arr)

            # Calculate Shannon Entropy and write to file
            entropy = round(self.shannon_entropy(filtered), 2)
            linetext = [col, str(entropy), str(n_gaps), str(self.n_seq-n_gaps), str(gap_fraction)]
            self.append_line_to_file(SE_outputfile, linetext)

            # Calculate similary using blosum62 matrix
            if self.blosum:
                seq = ''.join(filtered)
                blosum_score = self.calculate_similarity_score(seq)
                if np.isnan(blosum_score):
                    gap_fraction = np.nan
                linetext = [col, str(blosum_score), str(n_gaps), str(self.n_seq-n_gaps), str(gap_fraction)]
                self.append_line_to_file(blosum_outputfile, linetext)

        # If there are subfamilies provided; add to data frame in separate column
        if self.path_subfamilies:
            subgroups = extract_from_tsf(self.path_subfamilies, col_num=1, type=str)
            df['family'] = subgroups

            # Check that the minimum number of sequences per subfamily is reached
            subfam_counts = Counter(subgroups)
            for i in subfam_counts.items():
                if i[1] < self.threshold_sequences_per_subfamily:
                    sys.exit(f"Minimum number of sequences for {i[0]} is below thresholf of {self.threshold_sequences_per_subfamily}")

            # Calculate the Shannon Entropy per position in the MSA per subfamily and save in another file
            # 'shannon_entropy_{subgroupame}.txt'
            grouped = df.groupby('family')
            for name, group in grouped:

                print(f"processing {name}")

                # Create outputfile for SE
                SE_outputfile = os.path.join(self.outpath, f'shannon_entropy_{name}.txt')
                f = open(SE_outputfile, "w+")
                f.write('\t'.join(self.shannon_entropy_file_header) + '\n')
                f.close()

                if self.blosum:
                    blosum_outputfile = os.path.join(self.outpath, f'similarity_blosum62_{name}.txt')
                    f = open(blosum_outputfile, "w+")
                    f.write('\t'.join(self.similarity_file_header) + '\n')
                    f.close()

                for i in range(1, self.n_res + 1):
                    col = f'Pos_{i}'
                    arr = group[col].tolist()
                    n_gaps = arr.count(self.gap_symbol)  # Count and remove gaps
                    n_nogap = len(arr)
                    gap_fraction = round(n_gaps / len(arr), 2)
                    filtered = self.replace_gaps(arr)

                    # Calculate Shannon Entropy and write to file
                    entropy = round(self.shannon_entropy(filtered), 2)
                    linetext = [col, str(entropy), str(n_gaps), str(n_nogap), str(gap_fraction)]
                    self.append_line_to_file(SE_outputfile, linetext)

                    # Calculate similary using blosum62 matrix
                    if self.blosum:
                        seq = ''.join(filtered)
                        blosum_score = self.calculate_similarity_score(seq)
                        if np.isnan(blosum_score):
                            gap_fraction = np.nan
                        linetext = [col, str(blosum_score), str(n_gaps), str(self.n_seq-n_gaps), str(gap_fraction)]
                        self.append_line_to_file(blosum_outputfile, linetext)

        # Format results and store in separate file
        pos = self.positions_to_keep()
        res = self.extract_reference_residues() # In output store the residues of reference sequence
        self.format_results(pos, res)

        self.cleanup()  # Remove tab-sep file and non-formatted entropy results

        if self.calculate_pairwise_identity:  # Takes a long time to run
            # pairwise_seqidentity_matrix(self.path_msa)
            blosum62_scores_matrix(self.path_msa)

        print("Entropy calculatations finished")

    def replace_gaps(self, arr):
        """
        Replace gaps with a random character
        """
        output_list = []
        replace_char = 'X'
        for char in arr:
            if char == '-':
                output_list.append(replace_char)
                # Increment the character to use for the next '-' occurrence
                replace_char = chr(ord(replace_char) + 1)
            else:
                output_list.append(char)
        # print(output_list)
        return output_list

    def append_line_to_file(self, file, line):
        with open(file, 'a') as f:
            f.write('\t'.join(line) + '\n')

    def positions_to_keep(self):
        """
        Only keep positions that are ungapped in the reference file(s). Return list of these positions

        returns list of residue positions to keep (start from 1)
        """
        positionstokeep = set()
        # If reference ID(s) are provided, only take those positions, otherwise take all positions
        if self.reference_id:
            with open(self.path_msa) as handle: 
                for record in SeqIO.parse(handle, "fasta"):
                    if record.id in self.reference_id:
                        print("reference found: ", record.id)
                        for num, aa in enumerate(record.seq, 1):
                            if aa != self.gap_symbol:
                                positionstokeep.add(num)
            pos_list = list(positionstokeep)
            pos_list = [int(i) for i in pos_list]
            pos_list.sort()
        # Tab-separated file and no reference ID provided
        elif (not self.reference_id) and (self.msa_type == 'tab-separated'):
            with open(os.path.join(self.outpath, 'shannon_entropy.txt')) as f:
                lines = len(f.readlines()[1:])
                pos_list = [i for i in range(1, lines+1)]
        elif not (self.reference_id) and (self.msa_type == 'fasta'):
            with open(self.path_msa) as handle: 
                for record in SeqIO.parse(handle, "fasta"):
                    for num, aa in enumerate(record.seq, 1):
                        if aa != self.gap_symbol:
                            positionstokeep.add(num)
            pos_list = list(positionstokeep)
            pos_list = [int(i) for i in pos_list]
            pos_list.sort()
        return pos_list

    def extract_reference_residues(self):
        """
        d[ref_id] = [res1, res2, res3, res4]
        """
        result = {}
        
        # Check if there are multiple reference IDs (list) or just one (str)
        if isinstance(self.reference_id, list):
            for ref in self.reference_id:
                result[ref] = []
        elif isinstance(self.reference_id, str):
            result[self.reference_id] = []

        if not self.reference_id:
            return result

        with open(self.path_msa) as handle: 
                for record in SeqIO.parse(handle, "fasta"):
                    if record.id in self.reference_id:
                        count = 1
                        for aa in record.seq:
                            if aa == self.gap_symbol:
                                result[record.id].append(aa)
                            if not aa == self.gap_symbol:
                                result[record.id].append(aa + str(count))
                                count += 1 
        return result

    def format_results(self, positionstokeep, residuestokeep):
        for file in os.listdir(self.outpath):
            if 'shannon_entropy' or 'similarity_' in file:
                with open(os.path.join(self.outpath, file)) as f:
                    content = f.readlines()
                    outname_ext = os.path.splitext(file)
                    outname = outname_ext[0] + '_filtered.txt'
                    
                    with open(os.path.join(self.outpath, outname), 'w+') as out:
                        for num, line in enumerate(content):
                            if num in positionstokeep:

                                linelist = line.split()
                                alignment_pos = int(linelist[0].split('_')[1])

                                linelist.append(str(alignment_pos))  # Start at position 1 rather than 0
                                for k, v in residuestokeep.items():
                                    linelist.append(v[alignment_pos - 1])
                                out.write('\t'.join(linelist) + '\n')
                            elif num == 0:  # Header
                                linelist = line.split()
                                linelist.append('alignment_pos')
                                for k, v in residuestokeep.items():
                                    linelist.append(f"fullresidue_{k}")
                                out.write('\t'.join(linelist) + '\n')

    def make_df(self, M, ids):
        columns = [f'Pos_{i}' for i in range(1, self.n_res + 1)]
        df = pd.DataFrame(M, columns = columns)
        df['seqid'] = ids
        return df

    def format_data(self):        
        M = []
        ids = extract_from_tsf(self.path_msa_tsv, col_num=0, type=str)

        content = SE.read(self.path_msa_tsv)
        for line in content:
            line = line.split()
            aas = []
            for aa in line[1]:
                aas.append(aa)
            M.append(aas)

        SE.check_num_resulting_sequences(M)
        return M, ids

    def set_size_msa(self, length_sequences):
        self.n_res = length_sequences[0]
        self.n_seq = len(length_sequences)

    def cleanup(self):
        # Cleanup (remove tab-sep file and non-formatted entropy results)
        # - Tab-separated MSA file
        # - Non-filtered entropy values
        for file in os.listdir(self.outpath):
            if not '_filtered' in file:
                os.remove(os.path.join(self.outpath, file))
        # if not self.path_msa_tsv == self.path_msa:
        #     os.remove(self.path_msa_tsv)

    def check_size_msa(self):
        """
        Check whether length is same for all sequences in MSA.
        """
        lengths = extract_from_tsf(self.path_msa_tsv, col_num=1, type=str)
        lengths = [len(i) for i in lengths]
        check_size = SE.all_same(lengths)
        return check_size, lengths

    @staticmethod
    def get_tab_separated_fname(fp, extension="_tab-separated.txt"):
        return os.path.splitext(fp)[0] + extension

    @staticmethod
    def save_tab_separated_msa_totxt(fp):
        outname = SE.get_tab_separated_fname(fp)
        with open(outname, 'w+') as out:
            for record in SeqIO.parse(fp, "fasta"):
                out.write(record.id + '\t' + str(record.seq) + '\n')
    
    @staticmethod
    def all_same(items):
        return all(x == items[0] for x in items)

    @staticmethod
    def read(fname, headers=False):
        with open(fname, 'r') as f:
            if headers:
                content = f.readlines()[1:]
            else:
                content = f.readlines()
            return content

    @staticmethod
    def shannon_entropy(list_input):
        """
        Calculate Shannon's Entropy per column of an alignment (H=-\sum_{i=1}^{M} P_i\,log_2\,P_i)
        """

        unique_bases = set(list_input)
        entropy_list = []
        m = len(list_input)

        # Number of residues in column
        for base in unique_bases:
            n_i = list_input.count(base)  # Number of residues of type i
            p_i = n_i / float(m)  # n_i (Number of residues of type i) / M(Number of residues in column)
            entropy_i = p_i * (math.log(p_i, 2))
            entropy_list.append(entropy_i)
        sh_entropy = abs(-(sum(entropy_list)))
        return sh_entropy

    @staticmethod
    def check_num_resulting_sequences(M):
        if np.array(M).shape[0] == 0:
            sys.exit("No sequences in df")
        else:
            print(f"{str(np.array(M).shape[0])} sequences left")

    @staticmethod
    def is_fasta(filename):
        with open(filename, "r") as handle:
            fasta = SeqIO.parse(handle, "fasta")
            for i in fasta:
                if len(i.seq) < 1:
                    return False
                elif len(i.seq) > 1:
                    return True
            return any(fasta)  # False when `fasta` is empty, i.e. wasn't a FASTA file
        
    @staticmethod
    def natural_amino_acids():
        return 'ACDEFGHIKLMNPQRSTVWY'

    @staticmethod    
    def calculate_similarity_score(seq, method='BLOSUM62'):
        """
        Using BLOSUM62 matrix. Excluded diagonal and divide by length sequence
        """
        aligner = Align.PairwiseAligner()
        aligner.substitution_matrix = substitution_matrices.load(method)

        amino_acids = SE.natural_amino_acids()

        # Remove gaps from sequence
        # Remove non AA characters from sequence
        full_length = len(seq)
        seq = seq.replace('-', '')
        for char in seq:
            if char not in amino_acids:
                seq = seq.replace(char, '')
                
        num_gaps = full_length - len(seq)
        
        # Set gap penalty
        gap_penalty = 0.5 / full_length  # Penalty for each gap. 0.5 is max similarity score

        # Make matrix
        M = np.empty((len(seq), len(seq)), dtype='int')

        # Get upper/lower triangle (including diagonal) and calculate score
        for i in range(len(seq)):
            for num, j in enumerate(seq):
                score = aligner.substitution_matrix[str(seq[i]), str(seq[num])]
                M[i, num] = float(score)
      
        # Calculate score by summing values
        sum_tri = np.sum(M[np.triu_indices(len(seq), k=1)])  # Exclude diagonal
        sum_diagonal = np.trace(M)  
        if len(seq) == 1:  # Cannot (and should not) be calculated (Runtimerror)
            score = np.nan
        else:
            score = sum_tri / sum_diagonal  # Divide by diagonal for normalization between residues
            score = score - (gap_penalty * num_gaps)  # Subtract gap penalty
            score = round(score / full_length, 2)
        
        # If score is negative, set to 0
        if score < 0:
            score = 0

        return score
        

def main(args):
    
    instance = SE()
    instance.path_msa = args.msapath
    instance.outpath = args.outpath
    instance.path_subfamilies = args.subfamilypath
    instance.reference_id = args.reference_id
    
    if args.seqidentity:
        instance.calculate_pairwise_identity = True

    if args.blosum62:
        instance.blosum = True

    instance.run()

def read_arguments():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("-msa", "--msapath", required=True, help="file location of Multiple Sequence Alignment (FASTA format)")
    parser.add_argument("-sub", "--subfamilypath", required=True, help="file location of subfamily assignment (TXT file)")
    parser.add_argument("-out", "--outpath", required=True, help="output file location")
    parser.add_argument("-ref", "--reference_id", required=False, help="reference ID(s).", nargs="+")
    parser.add_argument("-id", "--seqidentity", help="Calculate pairwise sequence identity matrix", required=False)
    parser.add_argument("-b", "--blosum62", help="Calculate scores according to BLOSUM62 matrix", required=False, type=bool)
    parser.add_argument("-j", "--json", type=str, help="store command line argument as json. Provide folder where to store", required=False)
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    arguments = read_arguments()

    # Save command line argument in json dictionary
    if arguments.json:
        script = os.path.basename(__file__)
        save_commandline_args(script=script, args=vars(arguments), folder=arguments.json)

    main(arguments)