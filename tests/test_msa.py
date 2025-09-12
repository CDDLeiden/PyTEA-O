import pytest
import pathlib
import numpy as np
import pandas as pd

from src.utils.msa import MSA

@pytest.fixture
def test_file_path():
	return pathlib.Path(__file__).parent/'data'/'test.fasta.aligned'

@pytest.fixture
def temp_outdir():
	return 'temp'

def create_MSA_object():
	return MSA(msa_file=pathlib.Path(__file__).parent/'data'/'test.fasta.aligned',outdir='temp')

def test_MSA_creation(test_file_path,temp_outdir):
	
	msa_object = create_MSA_object()

	assert msa_object.msa_file == test_file_path
	assert msa_object.outdir == temp_outdir
	assert msa_object.num_sequences == 5
	assert msa_object.sequence_length == 6

def test_sequence_encoding():

	msa_object = create_MSA_object()

	msa = msa_object.msa

	actual_sequence = msa['seq3'].to_list()

	expected_sequence = [msa_object.char_to_int[x] for x in "ACDE-E"]

	assert actual_sequence == expected_sequence

def test_get_reference_indicies():

	msa_object = create_MSA_object()

	ref_idx = msa_object.get_reference_indices('seq3')

	assert np.array_equal(ref_idx,[0,1,2,3,5])

def test_get_residue_counts():

	msa_object = create_MSA_object()

	test_counts = {
		0: {i: 0 for i in range(24)},
		1: {i: 0 for i in range(24)},
		2: {i: 0 for i in range(24)},
		3: {i: 0 for i in range(24)},
		4: {i: 0 for i in range(24)},
		5: {i: 0 for i in range(24)},
	}

	# Manually filling in the correct counts
	test_counts[0][1] = 5
	
	test_counts[1][3] = 3
	test_counts[1][5] = 1
	test_counts[1][6] = 1
	
	test_counts[2][0] = 1
	test_counts[2][3] = 1
	test_counts[2][4] = 2
	test_counts[2][5] = 1
	
	test_counts[3][5] = 5
	
	test_counts[4][0] = 1
	test_counts[4][1] = 1
	test_counts[4][6] = 3
	
	test_counts[5][1] = 1
	test_counts[5][3] = 1
	test_counts[5][5] = 1
	test_counts[5][6] = 1
	test_counts[5][7] = 1

	test_counts = pd.DataFrame.from_dict(test_counts,orient='index')
	msa_counts = msa_object.get_residue_counts()

	test_index = test_counts.index.to_list()
	msa_index = msa_counts.index.to_list()

	assert np.array_equal(test_index,msa_index)

	test_columns = test_counts.columns.to_list()
	msa_columns = msa_counts.columns.to_list()

	assert np.array_equal(test_columns,msa_columns)
	
	identical_rows = 0
	for column in test_columns:
		letter = msa_object.int_to_char[column]
		test_data = test_counts[column].to_numpy()
		msa_data = msa_counts[column].to_numpy()
		if np.array_equal(test_data,msa_data): identical_rows += 1

	assert identical_rows == len(test_columns)

def test_calculate_sequence_difference():

	msa_object = create_MSA_object()

	expected_distances = np.array([
		[np.inf, 1.5, 2.0, 1.5, 2.0],
		[1.5, np.inf, 1.5, 1.0, 2.0],
		[2.0, 1.5, np.inf, 1.5, 1.5],
		[1.5, 1.0, 1.5, np.inf, 2.0],
		[2.0, 2.0, 1.5, 2.0, np.inf]
	])


	expected_distances = pd.DataFrame(expected_distances,columns=msa_object.accessions,index=msa_object.accessions)

	actual_distances = msa_object.calculate_sequence_difference(threads=2)

	assert np.array_equal(expected_distances.to_numpy(),actual_distances.to_numpy())