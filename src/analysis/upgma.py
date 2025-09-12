import pandas as pd
import numpy as np
import gc
import time
import typing
import argparse
import pathlib

from src.utils.msa import calculate_sequence_difference,load_msa
from src.utils.general import valid_file,valid_directory,parse_args,merge_commflags_with_kwargs

ARGUEMENTS = {
	'd':{
		'flag':'msa_dataframe',
		'required':True,
		'help':"Multiple Sequence Alignment dataframe created by load_dataframe",
		'type':pd.DataFrame,
		'CLI': False,
		'call': True
	},
	'o': {
		'flag': 'outdir',
		'help': "Output file location",
		'type': str,
		'default': "UPGMA_SUBFAMILY_GROUPINGS",
		'type': valid_directory,
		'CLI': True,
		'call': True
	},
	't': {
		'flag': 'threads',
		'help': "Number of threads",
		'type': int,
		'default': 1,
		'CLI': True,
		'call': True
	}
}

def UPGMA(distance_matrix:pd.DataFrame,outdir:pathlib.Path) -> dict:

	distance_matrix = distance_matrix + distance_matrix.T

	tree:dict = {0:distance_matrix.columns.to_list()}

	while(len(distance_matrix)>2):

		flattend_index =  np.argmin(distance_matrix)

		row_i = flattend_index//distance_matrix.shape[1]
		col_i = flattend_index%distance_matrix.shape[1]

		row = distance_matrix.columns[row_i]
		col = distance_matrix.columns[col_i]

		pairwise_dist = distance_matrix.loc[row,col]

		data_A = np.delete(distance_matrix.loc[row].to_numpy(),[row_i,col_i])
		weight_A = len(row.split(","))
		data_B = np.delete(distance_matrix.loc[col].to_numpy(),[row_i,col_i])
		weight_B = len(col.split(","))

		new_data = np.append((((data_A*weight_A) + (data_B*weight_B))/(weight_A+weight_B)),np.inf)

		distance_matrix = distance_matrix.drop(index=[row,col],columns=[row,col])

		joined_name = f"{row},{col}"

		tree[pairwise_dist/2] = distance_matrix.columns.to_list() + [joined_name]

		distance_matrix[joined_name] = new_data[:-1]
		distance_matrix.loc[joined_name] = new_data

		gc.collect()

	row,col = distance_matrix.columns.to_numpy()

	tree[distance_matrix.loc[row,col]] = [",".join(tree[0])]
	
	return tree

def run(args:argparse.Namespace=None,**kwargs) -> typing.Tuple[str,str]:

	args = merge_commflags_with_kwargs(
		cli_args=args,
		function_args=ARGUEMENTS,
		**kwargs
	)

	(outdir := args.outdir/'UPGMA').mkdir(mode=0o755,parents=True,exist_ok=True)

	distance_matrix:pd.DataFrame = calculate_sequence_difference(
		msa=args.msa_dataframe,
		threads=args.threads,
		outfile=f"{outdir}/upgma.dist"
	)

	tree = UPGMA(
		distance_matrix=distance_matrix,
		outdir=outdir
	)

	return tree

if __name__ == '__main__':

	run(args=parse_args(function_args=ARGUEMENTS))