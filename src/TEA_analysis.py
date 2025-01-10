#!/usr/bin/env python3

import os
import sys
import numpy as np
import pandas as pd
import plotly.express as px
import matplotlib.pyplot as plt
from sklearn import preprocessing
from utils import extract_from_tsf, read_alignment, load_pickle



class Analyse:
    """
    Set of functions to analyse Shannon Entropy results.
    """

    def __init__(self):
        self.mode = None
        self.se_folder = None
        self.reference = None

        self.gap_data = False  # Do not include gap data for plotting
        self.global_entropy_file = 'shannon_entropy.txt'
        self.summed_entropy_file = 'summed_subfamily_shannon_entropy.txt'
        self.average_entropy_file = 'average_shannon_entropy.txt'
        self.residue_counts_file = 'residue_counts.txt'

        self.data = {}

    def read_data(self):
        if self.reference:
            self.read_reference()
        if self.mode == "TEA":
            self.read_summed_entropy()
        elif self.mode == "TEAO":
            self.read_average_entropy()
        self.read_global_entropy()
        
        if not self.gap_data: # Remove gap positions
            gap_positions = [num for num, val in enumerate(self.data['residues']) if val == "-"]
            for index in sorted(gap_positions, reverse=True):
                for k, _ in self.data.items():
                    if isinstance(self.data[k], list) and len(self.data[k]) > 1:
                        try:
                            del self.data[k][index]
                        except IndexError:
                            print(f"Index {index} does not exist in {k}")


    def read_global_entropy(self):
        global_entropy = []
        msa_positions = []
        fraction_SE = []
        seq_at_pos = []
        gaps = []
        with open(os.path.join(self.se_folder, self.global_entropy_file), 'r') as f:
            content = f.readlines()[1:]  # Header
            for num, line in enumerate(content):
                line = line.strip().split()
                msa_pos = int(line[0])
                SE_glob = float(line[1])
                frac_SE = float(line[2])
                n_seq = int(line[3])
                gap = int(line[4])
                
                msa_positions.append(msa_pos)
                global_entropy.append(SE_glob)
                fraction_SE.append(frac_SE)
                seq_at_pos.append(n_seq)
                gaps.append(gap)

                if self.reference:
                    if not msa_pos in self.data['residue_MSA_positions']:
                        self.data['residues'].insert(num, "-")
                
        self.data['global_entropy'] = global_entropy
        self.data['msa_positions'] = msa_positions
        self.data['fraction_SE'] = fraction_SE
        self.data['N_seq_at_position'] = seq_at_pos
        self.data['N_gap_at_position'] = gaps


    def read_average_entropy(self):
        average_entropy = []
        with open(os.path.join(self.se_folder, self.average_entropy_file), 'r') as f:
            content = f.readlines()
            for line in content:
                line = line.strip().split()
                SE_avg = float(line[1])
                average_entropy.append(SE_avg)
        self.data['average_entropy'] = average_entropy


    def read_summed_entropy(self):
        summed_entropy = []
        with open(os.path.join(self.se_folder, self.summed_entropy_file), 'r') as f:
            content = f.readlines()
            for line in content:
                line = line.strip().split()
                SE_sum = float(line[1])
                summed_entropy.append(SE_sum)
        self.data['summed_entropy'] = summed_entropy


    def read_reference(self):
        residues = []
        residues_MSA_positions = []
        with open(os.path.join(self.se_folder, self.reference + '.shannon_entropy.txt'), 'r') as f:
            content = f.readlines()[1:] # Header
            for line in content:
                line = line.strip().split()
                pos = int(line[0]) + 1  # Start counting from 1
                msa_pos = int(line[1])
                res = line[2]
                residues.append(res + str(pos))
                residues_MSA_positions.append(msa_pos)
        self.data['residue_MSA_positions'] = residues_MSA_positions
        self.data['residues'] = residues


    def read_consensus_logo(self):
        pass


class Plot:
    """
    Create some of the plots for the TEA analysis
    """
    def __init__(self, data, mode):
        self.data = data
        self.mode = mode

    def global_vs_subfam(self, figsize=(5, 5)):
        f, ax = plt.subplots(figsize=figsize)
        if self.mode == "TEA":
            ax.scatter(self.data['summed_entropy'], self.data['global_entropy'])
        elif self.mode == "TEAO":
            ax.scatter(self.data['average_entropy'], self.data['global_entropy'])
        ax.set_xlabel('Average subfamily entropy')
        ax.set_ylabel('Global entropy')
        ax.legend()
        return f, ax
    
    def global_vs_subfam_interactive(self, width=800, height=800):
        if self.mode == "TEAO":
            fig = px.scatter(x=self.data['average_entropy'], 
                    y=self.data['global_entropy'],
                    hover_name=self.data['residues'],
                    width=width,
                    height=height)
            fig.update_layout(title='Average vs Global subfamily entropy', xaxis_title='Global entropy', yaxis_title='Average subfamily entropy')
        elif self.mode == "TEA":
            fig = px.scatter(x=self.data['summed_entropy'], 
                    y=self.data['global_entropy'],
                    hover_name=self.data['residues'],
                    width=width,
                    height=height)
            fig.update_layout(title='Summed vs Global subfamily entropy', xaxis_title='Global entropy', yaxis_title='Summed subfamily entropy')
        return fig