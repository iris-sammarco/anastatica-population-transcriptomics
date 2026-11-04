#!/bin/bash

#PBS -l walltime=4:00:00
#PBS -q default@pbs-m1.metacentrum.cz
#PBS -l select=1:ncpus=2:mem=40gb:scratch_local=20gb
#PBS -N 04_RDA_functional_enrichment_East_West
#PBS -j oe

# Usage:
# qsub 04_RDA_functional_enrichment_East_West_run.sh

set -euo pipefail

# Load conda/mamba environment
module add mambaforge

# Point to your conda environments directory
export CONDA_ENVS_PATH=~/.conda/envs/

# Activate the mamba environment where all the tools are installed
mamba activate dream_degs

# Point R temp to local scratch
export TMPDIR="${SCRATCHDIR}"

# Run the script
Rscript /Anastatica/scripts/04_RDA_functional_enrichment_East_West.R
