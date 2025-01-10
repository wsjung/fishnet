#!/bin/bash

######################################
# Master script for FISHNET pipeline #
######################################

### test parameters ###
trait="maleWC"
pvalFileDir=$( readlink -f ./test/${trait} )
pvalFileName=${pvalFileDir}/0-${trait}.csv
moduleFileDir=$( readlink -f ./test/ker_based )
numTests=$(( $(wc -l < "$pvalFileName") - 1 ))
geneColName="Genes"
pvalColName="p_vals"
bonferroni_alpha="0.05"
output_dir=$( readlink -f ./results/ )

### list of containers ###
# contains all python dependencies for fishnet
#   TODO: create single container with all python dependencies (include statsmodels)
#   TODO: add to biocontainers 
container_python="jungwooseok/dc_rp_genes:1.0"

################
#### PHASE 1 ###
################
## (1) nextflow original run
#echo "executing Nextflow MEA pipeline on original run"
#nextflow run ./scripts/phase1/nextflow/main.nf \
#    --trait $trait \
#    --moduleFileDir $moduleFileDir \
#    --numTests $numTests \
#    --pipeline $trait  \
#    --pvalFileName $pvalFileName \
#    --geneColName $geneColName \
#    --pvalColName $pvalColName  \
#    --bonferroni_alpha $bonferroni_alpha \
#    -c conf/fishnet.config
#
#
## (2) generate uniform p-values
#echo "generating uniformly distributed p-values"
#genes_filepath="/app/${pvalFileName#$(pwd)}"
#docker run -i --rm -v $(pwd):/app -u $(id -u):$(id -g) $container_python  /bin/bash -c \
#    "python3 /app/scripts/phase1/generate_uniform_pvals.py \
#       --genes_filepath $genes_filepath"
#echo "done"
#
#
## (3) nextflow random permutation run
#echo "executing Nextflow MEA pipeline on random permutations"
#pvalFileName=$( readlink -f ./test/${trait}RR/0-${trait}RR.csv )
#num_permutations="1"
#nextflow run ./scripts/phase1/nextflow/main.nf \
#    --trait "${trait}RR" \
#    --moduleFileDir $moduleFileDir \
#    --numTests $numTests \
#    --pipeline "${trait}RR"  \
#    --pvalFileName $pvalFileName \
#    --geneColName $geneColName \
#    --pvalColName $pvalColName  \
#    --bonferroni_alpha $bonferroni_alpha \
#    --random_permutation \
#    --numRP $num_permutations \
#    --GO_summaries_path "GO_summaries_RP" \
#    --masterSummaries_path "masterSummaries_RP" \
#    -c conf/fishnet.config \
#
## (4) compile results
## (4.1) original
#echo "compiling results"
#results_path="/app/${output_dir#$(pwd)}"
#summaries_path_original="${results_path}/masterSummaries/summaries/"
#docker run -i --rm -v $(pwd):/app -u $(id -u):$(id -g) $container_python /bin/bash -c \
#    "python3 /app/scripts/phase1/compile_results.py \
#        --dirPath $summaries_path_original \
#        --identifier $trait \
#        --output $results_path"
#
## (4.2) random permutation
#summaries_path_permutation="${results_path}/masterSummaries_RP/summaries/"
#docker run -i --rm -v $(pwd):/app -u $(id -u):$(id -g) $container_python /bin/bash -c \
#    "python3 /app/scripts/phase1/compile_results.py \
#        --dirPath $summaries_path_permutation \
#        --identifier ${trait}RR \
#        --output $results_path"
#echo "done"
#
## (5) organize phase 1 results
#echo "restructuring phase 1 results for input to phase 2"
## TODO: consider defining shell functions in separate file
## function to organize phase 1 nextflow pipeline results
#organize_nextflow_results() {
#    trait=$1
#    results_dir=$2
#    traitRR="${trait}RR"
#
#    # create outer trait directories
#    trait_dir="${results_dir}/${trait}/"
#    traitRR_dir="${results_dir}/${traitRR}/"
#    mkdir -p $trait_dir $traitRR_dir
#
#    # move over original results
#    mv "${results_dir}/GO_summaries" ${trait_dir}
#    mv "${results_dir}/masterSummaries" ${trait_dir}
#    mv "${results_dir}/master_summary_${trait}.csv" ${trait_dir}/master_summary.csv
#
#    # 3. move over permutation results
#    mv "${results_dir}/GO_summaries_RP" "${traitRR_dir}/GO_summaries"
#    mv "${results_dir}/masterSummaries_RP" "${traitRR_dir}/masterSummaries"
#    mv "${results_dir}/master_summary_${traitRR}.csv" ${traitRR_dir}/master_summary.csv
#}
#organize_nextflow_results $trait $output_dir
#echo "done"


###############
### PHASE 2 ###
###############
# (0) TODO: filter the summary file for original/permuted runs
#   - filter for 'True' in master_summary.csv -> master_summary_filtered.csv
#   - subset columns of master_summary_filtered.csv -> master_summary_filtered_parsed.csv
( head -n 1 "${output_dir}/${trait}/master_summary.csv"; grep 'True' "${output_dir}/${trait}/master_summary.csv" ) > "${output_dir}/${trait}/master_summary_filtered.csv"
cut -d ',' -f 1-8 "${output_dir}/${trait}/master_summary_filtered.csv" > "${output_dir}/${trait}/master_summary_filtered_parsed.csv"

( head -n 1 "${output_dir}/${trait}RR/master_summary.csv"; grep 'True' "${output_dir}/${trait}RR/master_summary.csv" ) > "${output_dir}/${trait}RR/master_summary_filtered.csv"
cut -d ',' -f 1-8 "${output_dir}/${trait}RR/master_summary_filtered.csv" > "${output_dir}/${trait}RR/master_summary_filtered_parsed.csv"

# (1) generate statistics for original run
genes_filedir="/app/${pvalFileDir#$(pwd)}"
module_filedir="/app/${moduleFileDir#$(pwd)}"
results_dir="/app/${output_dir#$(pwd)}"
echo "${results_dir}/${trait}/master_summary_filtered_parsed.csv"
for network in `ls ${moduleFileDir}`;
do
    docker run -i --rm -v $(pwd):/app -u $(id -u):$(id -g) $container_python /bin/bash -c \
        "python3 /app/scripts/phase2/dc_generate_or_statistics.py \
            --gene_set_path $genes_filedir \
            --master_summary_path ${results_dir}/${trait}/master_summary_filtered_parsed.csv \
            --trait $trait  \
            --module_path ${module_filedir}/${network} \
            --go_path ${results_dir}/${trait}/GO_summaries/${trait}/ \
            --study $trait \
            --output_path ${results_dir}/${trait}/results/raw/ \
            --network $network"
done

