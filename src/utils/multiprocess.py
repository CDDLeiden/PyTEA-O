import sys
import importlib

def get_multiprocessing_module():
	if 'darwin' in sys.platform:
		return importlib.import_module('multiprocess')
	return importlib.import_module('multiprocessing')