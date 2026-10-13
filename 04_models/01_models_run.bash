#!/bin/bash

#PBS -l walltime=12:00:00
#PBS -q default@pbs-m1.metacentrum.cz
#PBS -l select=1:ncpus=4:mem=20gb:scratch_local=20gb
#PBS -N 01_models
#PBS -j oe

# Usage:
# qsub 01_models_run.sh

# Load conda/mamba environment
module add mambaforge

# Point to your conda environments directory
export CONDA_ENVS_PATH= ~/.conda/envs/

# Activate the mamba environment where all the tools are installed
mamba activate dream_degs

# Point R temp to local scratch
export TMPDIR="${SCRATCHDIR}"

# Run the scripts
Rscript /Anastatica/scripts/01_single_models.r

Rscript /Anastatica/scripts/02_combined_models.r

Rscript /Anastatica/scripts/03_interaction_models.r


