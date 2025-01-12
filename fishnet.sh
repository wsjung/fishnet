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

###############
### PHASE 1 ###
###############
echo "
###############
### PHASE 1 ###
###############
"
# (1) nextflow original run
echo "# STEP 1.1: executing Nextflow MEA pipeline on original run"
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
echo "# STEP 1.2: generating uniformly distributed p-values"
genes_filepath="/app/${pvalFileName#$(pwd)}"
docker run -i --rm -v $(pwd):/app -u $(id -u):$(id -g) $container_python  /bin/bash -c \
    "python3 /app/scripts/phase1/generate_uniform_pvals.py \
       --genes_filepath $genes_filepath"
echo "done"


# (3) nextflow random permutation run
echo "# STEP 1.3: executing Nextflow MEA pipeline on random permutations"
echo "executing Nextflow MEA pipeline on random permutations"
pvalFileName=$( readlink -f ./test/${trait}RR/${trait}RR.csv )
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
echo "# STEP 4: compiling results"
# (4.1) original
echo "# STEP 4.1: original results"
results_path="/app/${output_dir#$(pwd)}"
summaries_path_original="${results_path}/masterSummaries/summaries/"
docker run -i --rm -v $(pwd):/app -u $(id -u):$(id -g) $container_python /bin/bash -c \
    "python3 /app/scripts/phase1/compile_results.py \
        --dirPath $summaries_path_original \
        --identifier $trait \
        --output $results_path"

# (4.2) random permutation
echo "# STEP 4.2: permutation results"
summaries_path_permutation="${results_path}/masterSummaries_RP/summaries/"
docker run -i --rm -v $(pwd):/app -u $(id -u):$(id -g) $container_python /bin/bash -c \
    "python3 /app/scripts/phase1/compile_results.py \
        --dirPath $summaries_path_permutation \
        --identifier ${trait}RR \
        --output $results_path"
echo "done"

# (5) organize phase 1 results
echo "# STEP 5: restructuring phase 1 results for input to phase 2"
# TODO: consider defining shell functions in separate file
# function to organize phase 1 nextflow pipeline results
organize_nextflow_results() {
    trait=$1
    results_dir=$2
    traitRR="${trait}RR"

    # create outer trait directories
    trait_dir="${results_dir}/${trait}/"
    traitRR_dir="${results_dir}/${traitRR}/"
    mkdir -p $trait_dir $traitRR_dir

    # move over original results
    mv "${results_dir}/GO_summaries" ${trait_dir}
    mv "${results_dir}/masterSummaries" ${trait_dir}
    mv "${results_dir}/master_summary_${trait}.csv" ${trait_dir}/master_summary.csv

    # 3. move over permutation results
    mv "${results_dir}/GO_summaries_RP" "${traitRR_dir}/GO_summaries"
    mv "${results_dir}/masterSummaries_RP" "${traitRR_dir}/masterSummaries"
    mv "${results_dir}/master_summary_${traitRR}.csv" ${traitRR_dir}/master_summary.csv
}
organize_nextflow_results $trait $output_dir
echo "done"


###############
### PHASE 2 ###
###############
echo "
###############
### PHASE 2 ###
###############
"
# (0) filter the summary file for original/permuted runs
echo "# STEP 2.0: filtering + parsing master_summary.csv files"
( head -n 1 "${output_dir}/${trait}/master_summary.csv"; grep 'True' "${output_dir}/${trait}/master_summary.csv" ) > "${output_dir}/${trait}/master_summary_filtered.csv"
cut -d ',' -f 1-8 "${output_dir}/${trait}/master_summary_filtered.csv" > "${output_dir}/${trait}/master_summary_filtered_parsed.csv"

( head -n 1 "${output_dir}/${trait}RR/master_summary.csv"; grep 'True' "${output_dir}/${trait}RR/master_summary.csv" ) > "${output_dir}/${trait}RR/master_summary_filtered.csv"
cut -d ',' -f 1-8 "${output_dir}/${trait}RR/master_summary_filtered.csv" > "${output_dir}/${trait}RR/master_summary_filtered_parsed.csv"

## (1) generate statistics for original run
echo "# STEP 2.1: generating statistics for original run"
genes_filedir="/app/${pvalFileDir#$(pwd)}"
module_filedir="/app/${moduleFileDir#$(pwd)}"
results_dir="/app/${output_dir#$(pwd)}"
for network in `ls ${moduleFileDir}/`;
do
    echo "Network: $network"
    network="${network%.*}" # remove extension
    docker run --rm -v $(pwd):/app -u $(id -u):$(id -g) $container_python /bin/bash -c \
        "python3 /app/scripts/phase2/dc_generate_or_statistics.py \
            --gene_set_path $genes_filedir \
            --master_summary_path ${results_dir}/${trait}/master_summary_filtered_parsed.csv \
            --trait $trait  \
            --module_path ${module_filedir}/${network}.txt \
            --go_path ${results_dir}/${trait}/GO_summaries/${trait}/ \
            --study $trait \
            --output_path ${results_dir}/${trait}/results/raw/ \
            --network $network"
done
echo "done"

# (2) generate statistics for permutation run
# TODO: rewrite in nextflow(?) for parallelism
# (2.1) TODO: prepare file with threshold:network pairs

# (2.2) generate statistics for permutation run
echo "# STEP 2.2: generating statistics for permutation runs"
genes_rpscores_filedir="/app/results/RPscores/${trait}RR/"
threshold_network_pairs=$( readlink -f ./test/slurm_thresholds_maleWC.txt )
while IFS=$'\t' read -r threshold network; do
    echo "Threshold $threshold, Network: $network"
    docker run --rm -v $(pwd):/app -u $(id -u):$(id -g) $container_python /bin/bash -c \
        "python3 /app/scripts/phase2/dc_generate_rp_statistics.py \
            --gene_set_path $genes_rpscores_filedir \
            --master_summary_path ${results_dir}/${trait}RR/master_summary_filtered_parsed.csv \
            --trait ${trait}RR \
            --module_path ${module_filedir}/${network}.txt \
            --go_path ${results_dir}/${trait}RR/GO_summaries/${trait}RR/ \
            --output_path ${results_dir}/${trait}RR/results/raw/ \
            --network $network \
            --threshold $threshold"
done < $threshold_network_pairs
echo "done"

# (3) summarize statistics
echo "# STEP 2.3: summarizing statistics"
for network in `ls ${moduleFileDir}/`;
do
    echo "Network: $network"
    network="${network%.*}" # remove extension
    docker run --rm -v $(pwd):/app -u $(id -u):$(id -g) $container_python /bin/bash -c \
        "python3 /app/scripts/phase2/dc_summary_statistics_rp.py \
            --trait $trait \
            --input_path $results_dir \
            --or_id $trait \
            --rr_id ${trait}RR \
            --input_file_rr_id ${trait}RR \
            --network $network \
            --output_path ${results_dir}/${trait}/summary/"
done

# (4) identify MEA passing genes
echo "# STEP 2.4: identify MEA-passing genes"
FDR_threshold=0.05 # TODO: add as input argument
percentile_threshold=0.99 # TODO: add as input argument
for network in `ls ${moduleFileDir}/`;
do
    echo "Network: $network"
    network="${network%.*}" # remove extension
    docker run --rm -v $(pwd):/app -u $(id -u):$(id -g) $container_python /bin/bash -c \
        "python3 /app/scripts/phase2/dc_identify_mea_passing_genes.py \
            --trait $trait \
            --geneset_input $genes_filedir \
            --FDR_threshold $FDR_threshold \
            --percentile_threshold $percentile_threshold \
            --network $network \
            --input_path ${results_dir}/${trait}"
done





echo "### FISHNET COMPLETE ###"
