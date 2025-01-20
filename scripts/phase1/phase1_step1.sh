#!/bin/bash

#SBATCH -J fishnet_phase1_step1
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH -o ./logs/fishnet_phase1_step1_%J.out

cd $1 # move to master script directory (important to set proper nextflow launchDir)

nextflow run ./scripts/phase1/nextflow/main.nf \
    --trait $TRAIT \
    --moduleFileDir $MODULEFILEDIR \
    --numTests $NUMTESTS \
    --pipeline $TRAIT  \
    --pvalFileName $PVALFILEPATH \
    --geneColName $GENECOLNAME \
    --pvalColName $PVALCOLNAME  \
    --bonferroni_alpha $BONFERRONI_ALPHA \
    -c $NXF_CONFIG
