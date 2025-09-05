import pandas as pd
import numpy as np

from src.utils.sequence import SequenceUtilities
from src.utils.msa import get_reference_indices

def write_references(reference_id:str,msa:pd.DataFrame,global_SE:pd.DataFrame,average_SE:np.ndarray,residues:np.ndarray,gaps:np.ndarray,outdir:str) -> None:

	print(f"\n\tSummarizing Two Entropy Analysis for reference sequence...")

	if not reference_id: 
		print(f"\n\t\t[N] Reference ID not provided, defaulting to first sequence in MSA [{msa.index[0]}]...")
		reference_id = msa.index[0]

	outfile = f"{outdir}/shannon_entropy.summary"

	reference_seq_to_msa_indexes = get_reference_indices(ref=reference_id,msa=msa)

	with open(outfile,'w') as OUT:
		OUT.write(f"## Sequence_Pos\tMSA_Pos\tResidue\tGlobal_Shannon_Entropy\tAverage_Shannon_Entropy\tSequences_at_pos\tNum_of_gaps\n")

		for res_index,msa_index in enumerate(reference_seq_to_msa_indexes):
			
			global_entropy = global_SE[0][msa_index]
			average_entropy = average_SE[0][msa_index]
			non_gapped_res = int(np.sum(residues[:,msa_index]))
			num_of_gaps = int(gaps[0][msa_index])

			## Residue index
			OUT.write(f"{res_index}")
			## MSA index
			OUT.write(f"\t{msa_index}")
			## Global entropy
			OUT.write(f"\t{global_entropy}")
			## Average entropy
			OUT.write(f"\t{average_entropy}")
			## Residue of reference
			OUT.write(f"\t{SequenceUtilities.AAs[np.where(msa.loc[reference_id,'array'].T[msa_index]==1)[0][0]]}")
			## Consenus
			residue_count = ";".join([f"{residue}:{int(count)}" for count,residue in sorted(list(zip(np.sum(msa.loc[msa.index,'array'],axis=0).T[msa_index],SequenceUtilities.AAs)),key=lambda x: x[0],reverse=True)])
			OUT.write(f"\t{residue_count}")
			## Non-gapped residues
			OUT.write(f"\t{non_gapped_res}")
			## Gapped residues
			OUT.write(f"\t{num_of_gaps}\n")

	return None

def write_shannon_entropy_temp_results(queue,temp_entropy_file) -> None:
	with open(temp_entropy_file,'a') as OUT:
		while True:
			data = queue.get()
			if not data:
				break
			key,entropy = data
			entropy = ",".join([f"{x}" for x in entropy[0]])
			OUT.write(f"{key}\t{entropy}\n")
	return None