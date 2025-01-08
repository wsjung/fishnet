#!/bin/bash

######################################
# Master script for FISHNET pipeline #
######################################

### test parameters ###
trait="maleWC"
pvalFileName=$( readlink -f ./test/${trait}/0-${trait}.csv )
moduleFileDir=$( readlink -f ./test/ker_based )
numTests=$(( $(wc -l < "$pvalFileName") - 1 ))
geneColName="Genes"
pvalColName="p_vals"
bonferroni_alpha="0.05"
output_dir=$( readlink -f ./results/ )

###############
### PHASE 1 ###
###############
# (1) nextflow original run
echo "executing Nextflow MEA pipeline on original run"
nextflow run ./scripts/phase1/nextflow/main.nf \
    --trait $trait \
    --moduleFileDir $moduleFileDir \
    --numTests $numTests \
    --pipeline $trait  \
    --pvalFileName $pvalFileName \
    --geneColName $geneColName \
    --pvalColName $pvalColName  \
    --bonferroni_alpha $bonferroni_alpha \
    -c conf/fishnet.config


# (2) generate uniform p-values
echo "generating uniformly distributed p-values"
genes_filepath="/app/${pvalFileName#$(pwd)}"
# contains all python dependencies for fishnet
# TODO: create single container with all python dependencies (include statsmodels)
# TODO: add to biocontainers 
container_python="jungwooseok/dc_rp_genes:1.0"
docker run -i --rm -v $(pwd):/app -u $(id -u):$(id -g) $container_python  /bin/bash -c \
    "python3 /app/scripts/phase1/generate_uniform_pvals.py \
       --genes_filepath $genes_filepath"
echo "done"


# (3) nextflow random permutation run
echo "executing Nextflow MEA pipeline on random permutations"
pvalFileName=$( readlink -f ./test/${trait}RR/0-${trait}RR.csv )
num_permutations="1"
nextflow run ./scripts/phase1/nextflow/main.nf \
    --trait "${trait}RR" \
    --moduleFileDir $moduleFileDir \
    --numTests $numTests \
    --pipeline "${trait}RR"  \
    --pvalFileName $pvalFileName \
    --geneColName $geneColName \
    --pvalColName $pvalColName  \
    --bonferroni_alpha $bonferroni_alpha \
    --random_permutation \
    --numRP $num_permutations \
    --GO_summaries_path "GO_summaries_RP" \
    --masterSummaries_path "masterSummaries_RP" \
    -c conf/fishnet.config \

# (4) compile results
# (4.1) original
echo "compiling results"
results_path="/app/${output_dir#$(pwd)}"
summaries_path_original="${results_path}/masterSummaries/summaries/"
docker run -i --rm -v $(pwd):/app -u $(id -u):$(id -g) $container_python /bin/bash -c \
    "python3 /app/scripts/phase1/compile_results.py \
        --dirPath $summaries_path_original \
        --identifier $trait \
        --output $results_path"

# (4.2) random permutation
summaries_path_permutation="${results_path}/masterSummaries_RP/summaries/"
docker run -i --rm -v $(pwd):/app -u $(id -u):$(id -g) $container_python /bin/bash -c \
    "python3 /app/scripts/phase1/compile_results.py \
        --dirPath $summaries_path_permutation\
        --identifier $trait \
        --output $results_path"
echo "done"





