# PyTEA-O

This repository contains a re-implementation of a Two-Entropy Analysis (TEA) in Python as described in [*Ye et al., 2008*](https://pubmed.ncbi.nlm.nih.gov/18304936/) for protein sequence variation analysis.

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

## Usage

For performing a TEA, the only input needed is a Multiple Sequence Alignment (MSA) in FASTA or CLUSTAL format.
To quickly get started, you can use the provided example MSA, this will allow you to verify that everything is working as expected.

First, you can create subgroups using the UPGMA algorithm 

```bash
python src/upgma_subfam_grouping.py \
    -m example/example.fasta \
    -o example/results
```

Next, the TEA calculations can be started using

```bash
python src/run_TEA.py \
    -m example/example.fasta \
    -o example/results \
    -r 0_0_0 \
    -f example/results/example_MSA.subfamiles
```

Where ```-r``` is the ID of the reference sequence and ```-f``` is the subfamily file that is created in the ```upgma_subfam_grouping.py```

## License

The software is licensed under the CC-BY-NC, which means it is free to use, remix, and build upon non-commercially as long as the copyright terms of the license are preserved. Commercial use is not permitted. You can view the LICENSE file for the full terms. If you have questions about the license or the use of the software in your organization, please contact Gerard J.P. van Westen at gerard@lacdr.leidenuniv.nl.
