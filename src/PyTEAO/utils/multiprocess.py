import sys
import importlib

import pathos

def Pool(threads:int=1):
	
	return pathos.multiprocessing.ProcessingPool(threads)