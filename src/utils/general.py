import os
import pathlib
import pandas as pd
import argparse
import inspect

def valid_file(path:str) -> pathlib.Path:
	p = pathlib.Path(path)
	if not p.is_file():
		raise argparse.ArgumentTypeError(f"File not found: {path}")
	return p

def parse_args(function_args:dict) -> argparse.Namespace:

	parser = argparse.ArgumentParser()

	for arg in function_args:

		meta:dict = function_args[arg]

		long_arg = meta.get('flag')
		required = meta.get('required')
		help = meta.get('help')
		type = meta.get('type')
		default = meta.get('default')
		choices = meta.get('choices')

		parser.add_argument(
			f"-{arg}",
			f"--{long_arg}",
			required=bool(required),
			help=help,
			type=type,
			default=default,
			choices=choices
		)

	return parser.parse_known_args()[0]

def merge_commflags_with_kwargs(cli_args:argparse.Namespace=None,function_args:dict|None=None,**kwargs):

	## Store default values for the function
	config = {} 
	for key in function_args.keys():
		if function_args[key].get('default'):
			config[function_args[key]['flag']] = function_args[key]['default']
		else:
			config[function_args[key]['flag']] = None
	
	## Update with command line arguments
	if cli_args:
		config.update(vars(cli_args))
	
	## Update with function call arguments
	config.update(kwargs)

	## Check for missing required arguments
	missing_args = [key for key in function_args.keys() if function_args[key].get('required') and function_args[key]['flag'] not in config]

	## Yell about missing args
	if missing_args:
		raise ValueError(f"Missing required keyword arguments: {', '.join(missing_args)}")
	
	## Make sure the args match their required types
	for spec in function_args.values():
	
		if not spec.get('type') or config[spec['flag']] is None:
			continue
		
		config[spec['flag']]  = spec['type'](config[spec['flag']]) 

	return argparse.Namespace(**config)

def get_file_name(file:str) -> str:
	"""
	Get the file name from the path without the extension
	"""
	return pathlib.Path(file).name.split('.')[0]

def create_directory(outdir:str) -> None:
	
	os.makedirs(outdir,0o750,exist_ok=True)

	return None

def read_subfamilies(subfamilies_file:str) -> dict:

	"""
	Reads subfamily file created by upgma_subfam_grouping.py

	#### Input
	*str* Filepath to subfamilies file


	#### Returns
	*dict* A hash of hashes containing all subfamilies for all branchpoints.
	The first set of keys are the branchpoints and the second is the group numbers for the subfamilies.
	The stored values are all the accessions for the corresponding group at the corresponding
	branchpoint.
	"""

	## Load in subfamily information; headers being keys and values being the subfamily it belongs to
	subfamilies = pd.read_csv(subfamilies_file,header=0,index_col=None,sep=';')

	families = {}

	## Each line represents a branch point; several sequences can be joined at the same branch point
	for index,branch in subfamilies.iterrows():

		families[index] = {}

		## Iterate over all sequences IDs, assigning them to their group numbers
		for key,value in branch.to_dict().items():

			## Initialize subfamliy at the current index
			if value not in families[index].keys():

				families[index][value] = []

			## Add the sequence ID to the proper subfamily
			families[index][value].append(key)

	return families