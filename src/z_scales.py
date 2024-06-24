import sys
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from prodec import ProteinDescriptors
from sklearn.preprocessing import MinMaxScaler


def calculate_protein_descriptor(aligned_sequences, descriptor):
    """
    Calculate protein descriptors of choice for unique sequences in the bioactivity dataset

    Parameters
    ----------
    aligned_sequences : list
        List of aligned sequences
    descriptor : str
        Descriptor to calculate

    Returns
    -------
    pandas.DataFrame
        Dataset with sequence and features for the protein descriptors of interest for sequences in the bioactivity dataset
    """
    desc_factory = ProteinDescriptors()
    prodec_descriptor = desc_factory.get_descriptor(descriptor)
    result = prodec_descriptor.pandas_get(aligned_sequences, gaps='omit')
    return result


def truncate_seq(msa, pos):
    """
    Truncate sequences in MSA to only include residues at positions of interest

    Parameters
    ----------
    msa : Bio.Align.MultipleSeqAlignment
        Multiple sequence alignment
    pos : list
        List of positions of interest

    Returns
    -------
    list
        List of truncated sequences
    """
    sequences = {}
    for record in msa:
        sequences[record.id] = []

    for p in pos:
        for record in msa:
            sequences[record.id].append(record.seq[p])

    truncated_seqs = []
    recids = []
    for k, v in sequences.items():
        sequences[k] = ''.join(v)
        truncated_seqs.append(sequences[k])
        recids.append(k)

    check_seq(truncated_seqs)  # Check if sequences contain only valid amino acids
    
    return truncated_seqs, recids


def check_seq(seq, print_invalid=False):
    """
    Check if sequence contains only valid amino acids

    Parameters
    ----------
    seq : list of sequences

    Returns
    -------
    list
    """
    aas = validaas()
    incorrect_positions = []
    for num, i in enumerate(seq):
        if not all(aa in aas for aa in i):
            if print_invalid:
                print('Invalid amino acid sequence found. Please check your input file.')
            incorrect_positions.append(num)
            # Find out which sequence is invalid
            for j in i:
                if j not in aas:
                    if print_invalid:
                        print('Invalid amino acid:', j)
                        print('Invalid sequence:', i)
                        print('Invalid sequence number:', num)
                        print('Please check your input file.')

    return incorrect_positions

def impute_seq(seqlist, type='remove'):
    """
    Impute missing residues in MSA

    Parameters
    ----------
    sequences : list of sequences
    pos : int
    type : str (default: 'remove') Remove gaps and replace with NaN value

    Returns
    -------
    list
        List of imputed sequences
    """
    # TODO FIX THIS
    columns = len(seqlist[0])
    colnames = [f'pos_{i}' for i in range(columns)]
    # TODO Add some statistics about what residues are imputed and what is the mode
    df_seqlist = pd.DataFrame(seqlist, columns=['seq'])
    df = df_seqlist['seq'].apply(lambda x: pd.Series(list(x)))
    df.replace('-', np.nan, inplace=True)
    return df
    # df['seq'] = df[:].agg(''.join, axis=1)
    # imputed_seqlist = df['seq'].tolist()
    # return imputed_seqlist


def scale_df(arr):
    df = pd.DataFrame(arr)
    df.head()
    scaler = MinMaxScaler()
    df_scaled = pd.DataFrame(scaler.fit_transform(df), columns=df.columns)
    return df_scaled


def std(df):
    return df.std(axis=0)


def reshape_df(df):
    arr = np.array(df)
    matrix = np.zeros((5, int(len(df)/5)))
    count = 0
    for y in range(int(len(df)/5)):
        for x in range(5):
            matrix[x][y] = arr[count]
            count += 1
    return matrix  # 5 different Z-scales (Sandberg)


def validaas():
    return ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'W', 'Y', 'V', 'M', 'N', 'P', 'Q', 'R', 'S', 'T']


def z_scales_heatmap(df, xlabels):
    fig, ax = plt.subplots(figsize=(20, 5))
    ylabels = ['Z1', 'Z2', 'Z3', 'Z4', 'Z5']  # Z-scales (Sandberg)
    sns.heatmap(df, ax=ax, cmap='Reds', annot=False)

    # put the major ticks at the middle of each cell
    ax.set_xticks(np.arange(df.shape[1]) + 0.5, minor=False)
    ax.set_yticks(np.arange(df.shape[0]) + 0.5, minor=False)

    # ax.set_xticks(np.arange(len(xlabels)))
    ax.set_xticklabels(xlabels, rotation=90)
    ax.set_yticklabels(ylabels, rotation=0)
    ax.xaxis.tick_top()
    ax.set_xlabel('Residue')
    ax.set_ylabel('Z-scale')
    ax.collections[0].colorbar.set_label('Standard deviation')  # Rename legend to std
    return fig, ax