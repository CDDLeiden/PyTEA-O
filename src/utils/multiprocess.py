import sys
import importlib
from src.io.output import write_shannon_entropy_temp_results

def get_multiprocessing_module():
	if 'darwin' in sys.platform:
		return importlib.import_module('multiprocess')
	return importlib.import_module('multiprocessing')