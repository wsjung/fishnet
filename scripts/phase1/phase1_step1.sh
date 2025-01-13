#!/bin/bash

#SBATCH -J fishnet_phase1_step1
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH -o ./logs/fishnet_phase1_step1_%J.out

echo "DEBUG: cd $1"
cd $1 # move to master script directory (important to set proper nextflow launchDir)

echo $(pwd)

nextflow run ./scripts/phase1/nextflow/main.nf \
    --trait $trait \
    --moduleFileDir $moduleFileDir \
    --numTests $numTests \
    --pipeline $trait  \
    --pvalFileName $pvalFileName \
    --geneColName $geneColName \
    --pvalColName $pvalColName  \
    --bonferroni_alpha $bonferroni_alpha \
    -c $NXF_CONFIG
