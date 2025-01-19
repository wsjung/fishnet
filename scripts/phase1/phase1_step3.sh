#!/bin/bash

#SBATCH -J fishnet_phase1_step3
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH -o ./logs/fishnet_phase1_step3_%J.out

cd $1 # move to master script directory (important to set proper nextflow launchDir)

nextflow run ./scripts/phase1/nextflow/main.nf \
    --trait $TRAITRR \
    --moduleFileDir $MODULEFILEDIR \
    --numTests $NUMTESTS \
    --pipeline $TRAITRR  \
    --pvalFileName $PVALFILEPATHRR \
    --geneColName $GENECOLNAME \
    --pvalColName $PVALCOLNAME  \
    --bonferroni_alpha $BONFERRONI_ALPHA \
    --random_permutation \
    --numRP $NUM_PERMUTATIONS \
    --GO_summaries_path "GO_summaries_RP" \
    --masterSummaries_path "masterSummaries_RP" \
    -c $NXF_CONFIG
