import os
import pathlib
import argparse
import hashlib
import tempfile
import pickle
import logging

logging = logging.getLogger(__name__)

def pickle_dump(obj:object,pickle_file:pathlib.Path):

	logging.info(f"Pickling {obj.__class__} object")

	pickle_file = pathlib.Path(pickle_file).resolve(strict=False)

	with tempfile.NamedTemporaryFile(dir=pickle_file.parent,delete=False) as tmp:
		pickle.dump(obj,tmp,protocol=pickle.HIGHEST_PROTOCOL)
		tmp_file = pathlib.Path(tmp.name)
	
	os.replace(tmp_file,pickle_file)

	logging.info(f"{obj.__class__} pickled succesfully.")


def pickle_load(pickle_file:pathlib.Path) -> object|None:

	pickle_file = valid_file(pickle_file,warn=False)

	if not pickle_file:
		return None

	with pickle_file.open("rb") as f:
		return pickle.load(f)
	

def get_file_hash(filepath:pathlib.Path):

	filepath = pathlib.Path(filepath).resolve(strict=False)

	with filepath.open("rb") as f:
		return hashlib.file_digest(f,"md5")


def valid_file(path:str,warn:bool=True) -> pathlib.Path|None:

	p = pathlib.Path(path).resolve(strict=False)
	
	if not p.is_file():
		
		if not warn:
			return None
		
		raise argparse.ArgumentTypeError(f"File not found: {path}")
	
	return p

def valid_directory(path:str) -> pathlib.Path:

	p = pathlib.Path(path).resolve(strict=False)

	try:
	
		## Make the directory and necessary parent directories
		p.mkdir(mode=0o755,parents=True,exist_ok=True)

		## Make sure directory is writable
		test_file = p/'.write_test'
		with open(test_file,'w') as f:
			f.write('test')

		## Remove test file
		test_file.unlink()

	except Exception as e:

		raise argparse.ArgumentTypeError(f"Directory {p} is not writable or cannot be created: {e}")

	return p

def binary(layout:str) -> str:

	from PyTEAO.utils.visualization import SUBPLOT_REGISTRY,import_all_subplots

	import_all_subplots()

	layout = layout[:len(SUBPLOT_REGISTRY.keys())]

	for i,x in enumerate(layout):

		if x not in "10":

			raise argparse.ArgumentTypeError(f"Invalid binary switch at position {i+1} ({"".join(["*" if i != j else y for j,y in enumerate(layout)])}), 0 or 1 is required.")

	return layout

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
		action = meta.get('action')

		kwargs = {
			'required':required is True,
			'help':help,
			'default':default,
		}

		if choices:
			kwargs['choices'] = choices

		if action:
			kwargs['action'] = action
		else:
			kwargs['type'] = type

		parser.add_argument(
			f"-{arg}",
			f"--{long_arg}",
			**kwargs
		)

	return parser.parse_known_args()[0]

def get_file_name(file:str) -> str:
	"""
	Get the file name from the path without the extension
	"""
	return pathlib.Path(file).name.split('.')[0]

def create_directory(outdir:str) -> None:
	
	os.makedirs(outdir,0o750,exist_ok=True)

	return None