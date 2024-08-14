#!/usr/bin/env python3

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn import preprocessing
from utils import extract_from_tsf, read_alignment, load_pickle


class TEA_analysis:
    """
    Set of functions to analyse Shannon Entropy results.

    """
    def __init__(self):
        
        # File paths
        self.input_directory = None
        self.global_entropy_file = 'shannon_entropy_filtered.txt'
        self.global_similarity_file = 'similarity_blosum62_filtered.txt'

        self.domains_file = None
        self.annotation_file = None
        self.output_location = None
        self.subgroup_seqid_mapping = None  # Dictionary of format d[prot] = [seqid1, seqid2 etc ...]
        self.scale = True # Scale the entropy/similarity values to 0-1

        # Calculated from input
        self.domains = None
        self.annotation = None
        self.global_entropy = dict()
        self.subfam_entropy = dict()
        self.proteins = None  # The protein names are the names provided in the subrgoups file of the shannon_entropy.py script
        self.similarity = False

    def run(self):
        
        # Add global entropy to dict
        self.read_global_entropy()

        # Loop over subfam entropy files and add to entropy dictionary
        entropy_subfam_files = [file for file in os.listdir(self.input_directory) if (file != self.global_entropy_file) and 'shannon_entropy_' in file]  # Subfam files are saved as shannon_entropy_dnae1_filtered.txt > dnae1 is extracted
        similarity_subfam_files = [file for file in os.listdir(self.input_directory) if (file != self.global_similarity_file) and 'similarity_blosum62_' in file]
        if len(similarity_subfam_files) >= 1:
            print("Similarity files found, adding similarity to entropy dictionary")
            self.similarity = True

        self.proteins = [i.split('_')[2] for i in entropy_subfam_files]

        self.read_subfam_entropy(entropy_subfam_files)
        
        # Read similarity files
        if self.similarity:
            self.read_subfam_similarity(similarity_subfam_files)
        
        self.add_subfam_entropies()

        if self.domains_file:
            self.read_domains()

        if self.annotation_file:
            self.read_annotation()

        if self.scale:
            self.scale_results()

    def make_df(self):
        df = pd.DataFrame()
        for key, value in self.global_entropy.items():
            if not len(value) == 0:
                df[key] = value
        for key, value in self.subfam_entropy.items():
            if not len(value) == 0:
                df[key] = value
        return df
    
    def scale_results(self):
        scaler = preprocessing.MinMaxScaler(feature_range=(0, 1))
        for k, v in self.subfam_entropy.items():
            if 'entropy' in k:
                tmp = scaler.fit_transform(np.array(self.subfam_entropy[k]).reshape(-1, 1)).tolist()
                tmp2 = [item for sublist in tmp for item in sublist]
                self.subfam_entropy[k] = tmp2
        for k, v in self.global_entropy.items():
            if 'entropy' in k:
                tmp = scaler.fit_transform(np.array(self.global_entropy[k]).reshape(-1, 1)).tolist()
                tmp2 = [item for sublist in tmp for item in sublist]
                self.global_entropy[k] = tmp2

    def add_subfam_entropies(self):
        added_subfam_entropy = []
        added_subfam_similarity = []
        count = 0
        for prot in self.proteins:
            added_subfam_entropy.append(self.subfam_entropy[f"{prot}_entropy"])
            if not self.subfam_entropy.get(f"{prot}_similarity") == None:
                added_subfam_similarity.append(self.subfam_entropy[f"{prot}_similarity"])
            count += 1
        summed_subfam_SE = [sum(x) for x in zip(*added_subfam_entropy)]
        summed_subfam_similarity = [sum(x) for x in zip(*added_subfam_similarity)]
        average_subfam_SE = [i/count for i in summed_subfam_SE]
        average_subfam_similarity = [i/count for i in summed_subfam_similarity]
        self.subfam_entropy['summed_subfam_entropy'] = summed_subfam_SE
        self.subfam_entropy['average_subfam_entropy'] = average_subfam_SE
        self.subfam_entropy['summed_subfam_similarity'] = summed_subfam_similarity
        self.subfam_entropy['average_subfam_similarity'] = average_subfam_similarity

    def find_highly_gapped_positions(self, max_fraction=0.8):
        # Remove highly filtered positions from alignment
        # Global entropy and subfam entropy
        filtered_residues = []
        for prot in self.proteins:
            for num, val in enumerate(self.subfam_entropy[f"{prot}_gap_fraction"]):
                if (val > max_fraction) or (val == np.nan):
                    filtered_residues.append(num)
        for num, val in enumerate(self.global_entropy["gap_fraction"]):
            if (val > max_fraction) or (val == np.nan):
                        filtered_residues.append(num)
        filtered_residues = list(set(filtered_residues))
        return filtered_residues  # Residues to remove!

    def filter_residues(self, indexes):
        for index in sorted(indexes, reverse=True):
            for k, _ in self.global_entropy.items():
                # Check if the key is a list (reference sequence is not a list)
                if isinstance(self.global_entropy[k], list) and len(self.global_entropy[k]) > 1:
                    try:
                        del self.global_entropy[k][index]
                    except IndexError:
                        print(f"Index {index} does not exist in {k}")
            for k, _ in self.subfam_entropy.items():
                if isinstance(self.subfam_entropy[k], list) and len(self.subfam_entropy[k]) > 1:
                    try:
                        del self.subfam_entropy[k][index]
                    except IndexError:
                        print(f"Index {index} does not exist in {k}")

    def mpl_TEA_plot_2D(self, x, y, xlabel, ylabel):
        f, ax = plt.subplots(figsize=(5, 5))
        ax.scatter(x=x, y=y)
        ax.set_aspect('equal', adjustable='box')
        ax.set_ylabel(ylabel)
        ax.set_xlabel(xlabel)
        return f, ax

    # def plot_alignment_entropy(self, alpha=0.5):
    #     f, axs = plt.subplots(len(self.proteins) + 1)
    #     f.suptitle('Vertically stacked subplots')
    #     axs[0].plot(self.global_entropy['residues'], self.global_entropy['global_entropy'], label = 'global entropy', alpha=alpha)
    #     for num, prot in enumerate(self.proteins):
    #         axs[num + 1].plot(self.global_entropy['residues'], self.subfam_entropy[f'{prot}_entropy'], label = f"{prot}", alpha=alpha)
    #     axs.legend()
    #     axs.set_xlabel('alignment position')
    #     axs.set_ylabel('SE')
    #     return f, axs

    def plot_alignment_entropy(self, alpha=0.5):
        f, ax = plt.subplots(figsize=(15, 5))
        ax.plot(self.global_entropy['residues'], self.global_entropy['global_entropy'], label = 'global entropy', alpha=alpha)
        for num, prot in enumerate(self.proteins):
            ax.plot(self.global_entropy['residues'], self.subfam_entropy[f'{prot}_entropy'], label = f"{prot}", alpha=alpha)
        ax.legend()
        ax.set_xlabel('alignment position')
        ax.set_ylabel('SE')
        return f, ax

    def read_subfam_entropy(self, subfam_files):
        for prot, file in zip(self.proteins, subfam_files):
            self.subfam_entropy[f"{prot}_entropy"] = extract_from_tsf(os.path.join(self.input_directory, file), 1, type=float, begin=1)
            self.subfam_entropy[f"{prot}_n_gaps"] = extract_from_tsf(os.path.join(self.input_directory, file), 2, type=int, begin=1)
            self.subfam_entropy[f"{prot}_gap_fraction"] = extract_from_tsf(os.path.join(self.input_directory, file), 4, type=float, begin=1)
            self.subfam_entropy[f"{prot}_residues"] = extract_from_tsf(os.path.join(self.input_directory, file), 5, type=int, begin=1)

    def read_subfam_similarity(self, subfam_files):
            for prot, file in zip(self.proteins, subfam_files):
                self.subfam_entropy[f"{prot}_similarity"] = extract_from_tsf(os.path.join(self.input_directory, file), 1, type=float, begin=1)
            
    def read_global_entropy(self):
        self.global_entropy['global_entropy'] = extract_from_tsf(os.path.join(self.input_directory, self.global_entropy_file), 1, type=float, begin=1)        
        self.global_entropy['n_gaps'] = extract_from_tsf(os.path.join(self.input_directory, self.global_entropy_file), 2, type=int, begin=1)
        self.global_entropy['gap_fraction'] = extract_from_tsf(os.path.join(self.input_directory, self.global_entropy_file), 4, type=float, begin=1)
        self.global_entropy['residues'] = extract_from_tsf(os.path.join(self.input_directory, self.global_entropy_file), 5, type=int, begin=1)

        # Load global similarity if it exists
        if os.path.isfile(os.path.join(self.input_directory, self.global_similarity_file)):
            self.global_entropy['global_similarity'] =  extract_from_tsf(os.path.join(self.input_directory, self.global_similarity_file), 1, type=float, begin=1)

        # TODO MOVE TO SEPARATE FUNCTION
        with open(os.path.join(self.input_directory, self.global_entropy_file), 'r') as f:
            content = f.readlines()
            header = content[0]
            header = header.split()
        for num, ref in enumerate(header[6:]):
            refidlist = ref.split('_')[1:]
            refid ='_'.join(refidlist)
            self.global_entropy[f"{refid}"] = extract_from_tsf(os.path.join(self.input_directory, self.global_entropy_file), num + 6, type=str, begin=1)

    def read_domains(self):
        """
        Domainsfile in format:
        
        ';' separated CSV file
        for every residue write resnum;restype;domain such as:
        
            1;M;no domain found 
            1000;W;domain1 

        header should be the following:

            position; residue; domain
        
        """
        df = pd.read_csv(self.domains_file, comment='#', sep=';')
        fullres_domains = {}  # d[M1] = ['PHP']
        for _, row in df.iterrows():
            resnum = int(row['position'])
            aa = row[' residue']
            aa = aa.split()[0]
            fullres_domains[aa + str(resnum)] = row[' domain']
        self.domains = fullres_domains
        return fullres_domains

    def read_annotation(self):
        """
        Annotationfile in format:
        
        ';' separated CSV file
        for every residue writeres_num;aa;annotation such as:
        
            860;G;Resistance-causing mutation in S. aureus

        write 1 line per resnum for all residue numbers. 

        header should be the following:

            res_num;aa;residue_function
        """
        df = pd.read_csv(self.annotation_file, sep='[;,]', error_bad_lines=False)
        fullres_annotation = {}  # d[M1] = ['annotation']
        for index, row in df.iterrows():
            resnum = int(row['res_num'])
            aa = row['aa']
            aa = aa.split()[0]
            if not pd.isna(row['residue_function']):
                fullres_annotation[aa + str(resnum)] = row['residue_function']
        self.annotation = fullres_annotation
        return fullres_annotation

    @staticmethod
    def extract_residues_from_position(msa_path: str, alignment_pos: int, mapping_fp: str):
        # Alignment pos starting from 0
        # Read alignment
        # Get all residue from msa at that position
        residues = {}
        mapping = load_pickle(mapping_fp)  # Mapping between alignment and subgroup, e.g. ID1 prot1, ID2, prot2
        alignment = read_alignment(msa_path)
        for k, v in mapping.items():
            residues[k] = []

        for record in alignment:
            for k, v in mapping.items():
                if record.id in v:
                    residues[k].append(record.seq[alignment_pos])  # Count in file starting from 1
        # d[prot1] = [A, A, A]
        # d[prot2] = [K, K, K]       
        return residues

    @staticmethod
    def plot_residues_from_position(msa_path: str, alignment_pos: int, mapping_fp: str):
        residues = TEA_analysis.extract_residues_from_position(msa_path, alignment_pos, mapping_fp)
        # TODO Add grouped histogram
        pass

    @staticmethod
    def mpl_hist(x):  # Matplotlib
        f, ax = plt.subplots(figsize=(5, 2))  # figsize=(10, 10)
        ax.hist(x, density=False, bins=60)  # density=False would make counts
        ax.set_ylabel('Count')
        ax.set_xlabel('Shannon Entropy')
        return f, ax
