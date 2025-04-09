#!/bin/bash

SCRIPT_DIR="$(dirname "$0")"
codebasedir="${SCRIPT_DIR}/../src"
configfile="${SCRIPT_DIR}/../config/graph_results.cfg"
msafile="${SCRIPT_DIR}/example_MSA.clustal"
outputdir="${SCRIPT_DIR}/results"
refid="0_0_0"

# Create output directory
mkdir -p ${outputdir}

# Create one-liner alignment
echo "Create one-liner .clustal alignment"
msafolder=$(dirname "${msafile}")
perl ${codebasedir}/make_aln_oneliner.pl \
    -i ${msafile} \
    -o ${msafolder}

# Create subgroups
echo "Make subgroups"
python ${codebasedir}/upgma_subfam_grouping.py \
    -m ${msafile}

# Run TEA-O calculations
echo "Run TEA calculations"
python ${codebasedir}/shannon_entropy.py \
    -m "example_MSA/example_MSA.olaln" \
    -o ${outputdir} \
    -p "TEAO" \
    -r ${refid} \
    -s "example_MSA_grouping.txt" \
    -t 12

# Calculate z-scales
echo "Calculate z-scales"
python ${codebasedir}/z_scales.py \
    -m "example_MSA/example_MSA.olaln" \
    -o ${outputdir}

# plot the results
echo "Plot the results"
python ${codebasedir}/graph_results.py \
    -s ${outputdir}/${refid}.shannon_entropy.txt \
    -c ${outputdir}/consensus_logo.txt \
    -a ${outputdir}/teao_average_shannon_entropy.txt \
    -z ${outputdir}/zscales.txt \
    -o ${outputdir} \
    -l ${configfile}