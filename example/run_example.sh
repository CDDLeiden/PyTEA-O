#!/bin/bash

SCRIPT_DIR="$(dirname "$0")"
codebasedir="${SCRIPT_DIR}/../src"
configfile="${SCRIPT_DIR}/../config/graph_results.cfg"
msafile="${SCRIPT_DIR}/example.fasta"
outputdir="${SCRIPT_DIR}/results"
refid="0_0_0"

# Create subgroups
echo "Make subgroups"
python ${codebasedir}/upgma_subfam_grouping.py \
    -m ${msafile} \
    -o ${outputdir}

# Run TEA-O calculations
echo "Run TEA calculations"
python ${codebasedir}/run_TEA.py \
    -m ${msafile} \
    -o ${outputdir} \
    -r ${refid} \
    -t 12 \
    -f ${outputdir}/example_MSA.subfamilies