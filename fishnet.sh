#!/bin/bash

######################################
# Master script for FISHNET pipeline #
######################################

### test parameters ###
trait="maleWC"
pvalFileName=test/maleWC/0-maleWC.csv
moduleFileDir=test/ker_based
numTests=$(( $(wc -l < "$pvalFileName") - 1 ))
geneColName="Genes"
pvalColName="p_vals"

nextflow run ./scripts/phase1/nextflow/main.nf \
    --trait $trait \
    --moduleFileDir $moduleFileDir \
    --numTests $numTests \
    --pipeline $trait  \
    --pvalFileName $pvalFileName \
    --geneColName $geneColName \
    --pvalColName $pvalColName  \
    -c conf/fishnet.config
