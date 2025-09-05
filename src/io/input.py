

def read_temp_tea_file(temp_file:str) -> dict:

	results = {}

	with open(temp_file,'r') as IN:
		for line in IN:
			line = line.strip()
			if line == "":
				continue
			key,data = line.split("\t")
			results[int(key)] = [float(x) for x in data.split(",")]

	return results