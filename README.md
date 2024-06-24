# TEA

This repository contains an implementation of the Two-Entropy Analysis (TEA) in Python as described in [*Ye et al., 2008*](https://pubmed.ncbi.nlm.nih.gov/18304936/)

## Setup

### Prerequisites
- [Anaconda](https://www.anaconda.com/products/distribution) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html) installed on your system.

### Creating a Conda Environment

1. **Clone the repository:**
    ```sh
    git clone https://gitlab.services.universiteitleiden.nl/cdd/tea.git
    cd tea
    ```

2. **Create a new Conda environment:**
    ```sh
    conda create --name TEA
    ```

3. **Activate the Conda environment:**
    ```sh
    conda activate TEA
    ```

4. **Install the required dependencies:**
    ```sh
    conda install --file environment.yaml
    ```

5. **(Optional) Deactivate the Conda environment:**
    ```sh
    conda deactivate
    ```

## Quick Start

For performing a TEA, you will need a MSA in fasta format.

Here we are using two test directories, each containing a Multiple Sequence Alignment of different bacterial DNA polymerases. In the steps below we describe how to run the scripts for a protein-based as well as lineage-based analysis.

*protein-based*

Here, we are going to analyse two different DNA polymerases, DnaE1 and DnaE2. In the example a MSA was made containing sequences of both proteins. the files dnae1.fasta and dnae2.fasta contain the sequences that belong to each of the two polymerases

*lineage-based*

Here, we are analysing two polymerases (dnae and polC), based on their lineage. The file dnae_polc_alignment.afa contains the MSA in fasta format



### 1. Obtain the lineage

For every sequence in the MSA obtain the lineage using the Entrez API. Running this script takes some time as for every sequence the Entrez API will be called. The output will be two files with [MSA-filename]_filtered.fasta containing the sequences for which the lineage could be found and [MSA-filename]_lineage.csv which contains the lineage information. These files can be used for dividing the sequences in subgroups

```sh
python src/obtain_lineage.py \
    -i test-files/lineage-based/dnae_polc_alignment.afa \
    -o test-files/lineage-based \
    -e Your.Name.Here@example.org
```

### 2. Divide the sequences into subgroups (lineage-based or protein-based)

The sequences will be divided into subgroups, based on the lineage with the script '1a_subgroups_lineage.py' or protein-based using '1b_subgroups_protein-based.py'

```sh
python src/1a_subgroups_lineage.py \
    -l test-files/lineage-based/dnae_polc_alignment_lineage.csv \
    -s test-files/lineage-based \
    -m test-files/lineage-based/dnae_polc_alignment_filtered.fasta
    -r WP_161869988.1
```

```sh
python src/1b_subgroups_protein-based.py \
    -i test-files/protein-based/dnae_polc_alignment_filtered.fasta \
    -o test-files/protein-based \
    -s test-files/protein-based/dnae1.fasta test-files/protein-based/dnae2.fasta
```

### 3. Running the entropy calculations

_protein-based_

```sh
python src/shannon_entropy.py \
    -msa test-files/protein-based/dnae_polc_alignment_filtered.fasta \
    -sub test-files/protein-based/subgroups.txt \
    -out test-files/protein-based/subgroups.txt \
    -ref P9WNT7
```

_lineage-based_

Here, wer are going to calculate the entropy per family

```sh
python src/shannon_entropy.py \
    -msa test-files/lineage-based/family_dnae_polc_alignment_filtered.fasta \
    -sub test-files/lineage-based/species-based/subfamily_family.txt \
    -out test-files/lineage-based/species-based/results \
    -ref WP_161869988.1
```

Apart from the Shannon Entropy, you can also calculate the 'Shannon Similarity', making use of the BLOSUM matrix.

### 4. Visualize the results

Results are visualized in the jupyther notebooks:

- test-analysis\test-analysis-lineage-based.ipynb
- test-analysis\test-analysis-protein-based.ipynb


## License

The software is licensed under the CC-BY-NC, which means it is free to use, remix, and build upon non-commercially as long as the copyright terms of the license are preserved. Commercial use is not permitted. You can view the LICENSE file for the full terms. If you have questions about the license or the use of the software in your organization, please contact Gerard J.P. van Westen at gerard@lacdr.leidenuniv.nl.

