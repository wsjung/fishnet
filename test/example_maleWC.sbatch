#!/bin/bash

#SBATCH -J fishnet_example_maleWC
#SBATCH --mem=4G
#SBATCH -o ./logs/fishnet_example_maleWC_%J.out

# load required modules
eval $(spack load --sh singularityce@3.8.0)
eval $(spack load --sh nextflow@22.10.4)

# recommended environmental variables
export SINGULARITY_CACHEDIR=""
export NXF_CONDA_CACHEDIR=""
export TMPDIR=""
export NXF_TEMP=""

study="./data/pvals/exampleOR/"
study_random="./data/pvals/exampleRR/"
modules="./data/modules/ker_based/"
num_permutations=500

./fishnet_multi.sh \
    --study $study \
    --study-random $study_random \
    --modules $modules \
    --singularity \
    --conda \
    --conda_env fishnet \
    --num-permutations $num_permutations

