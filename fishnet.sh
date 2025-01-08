#!/bin/bash

######################################
# Master script for FISHNET pipeline #
######################################

### test parameters ###
trait="maleWC"
pvalFileName=$( readlink -f ./test/maleWC/0-maleWC.csv )
moduleFileDir=$( readlink -f ./test/ker_based )
numTests=$(( $(wc -l < "$pvalFileName") - 1 ))
geneColName="Genes"
pvalColName="p_vals"
bonferroni_alpha="0.05"

###############
### PHASE 1 ###
###############
# (1) nextflow original run
nextflow run ./scripts/phase1/nextflow/main.nf \
    --trait $trait \
    --moduleFileDir $moduleFileDir \
    --numTests $numTests \
    --pipeline $trait  \
    --pvalFileName $pvalFileName \
    --geneColName $geneColName \
    --pvalColName $pvalColName  \
    --bonferroni_alpha $bonferroni_alpha \
    -c conf/fishnet.config \

# (2) generate uniform p-values
genes_filepath="/app/${pvalFileName#$(pwd)}"
container_phase1_step2="jungwooseok/dc_rp_genes:1.0"    # TODO: add to biocontainers
docker run -i --rm -v $(pwd):/app -u $(id -u):$(id -g) jungwooseok/dc_rp_genes:1.0 /bin/bash -c \
    "python3 /app/scripts/phase1/generate_uniform_pvals.py \
       --genes_filepath $genes_filepath"
