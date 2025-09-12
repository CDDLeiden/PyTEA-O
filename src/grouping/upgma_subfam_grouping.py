#!/usr/bin/env python3

from os.path import isfile
import numpy as np
import pandas as pd
from warnings import simplefilter
import pathlib

simplefilter(action='ignore',category=pd.errors.PerformanceWarning)

def generate_subgroups(tree:dict,outfile:pathlib.Path)-> None:
	
	subfamilies = {}
	for branchpoint, distance in enumerate(sorted(tree.keys())):
		subfamilies[branchpoint] = {}
		for group_num,group in enumerate(sorted(tree[distance])):
			for member in group.split(","):
				subfamilies[branchpoint][member] = hex(group_num)

	with outfile.open("w") as SUBFAM:
		headers = ";".join(tree[0])
		SUBFAM.write(f"{headers}\n")
		for branch in subfamilies:
			grouping_line = ";".join([subfamilies[branch][member] for member in sorted(subfamilies[branch].keys())])
			SUBFAM.write(f"{grouping_line}\n")

	return None