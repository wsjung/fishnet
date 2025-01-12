#!/bin/bash -ue

######################################
# Master script for FISHNET pipeline #
######################################

usage() {
cat << EOF
Master script for the FISHNET pipeline
Usage: fishnet.sh [options]

  Options:
    -h, --help
      Prints the command usage
      Default: false
    --test
      Runs the test configuration
      Default: true
    --skip-stage-1
        Skips stage 1 of FISHNET
        Default: false
    --skip-stage-2
        Skips stage 2 of FISHNET
        Default: false
    --thresholding-alternative
        Configures stage 2 to run alternative thresholding mechanism
        Default: false (runs default thresholding mechanism)
EOF
}

TEST_MODE=false
SKIP_STAGE_1=false
SKIP_STAGE_2=false
THRESHOLDING_MODE_DEFAULT="default"
THRESHOLDING_MODE_ALTERNATIVE="alternative"
THRESHOLDING_MODE=$THRESHOLDING_MODE_DEFAULT

while [[ $# -gt 0 ]]; do
    case "$1" in
        -h|--help)
            usage
            exit 0
            ;;
        --test)
            TEST_MODE=true
            shift
            ;;
        --skip-stage-1)
            SKIP_STAGE_1=true
            shift
            ;;
        --skip-stage-2)
            SKIP_STAGE_2=true
            shift
            ;;
        --thresholding-alternative)
            THRESHOLDING_MODE=$THRESHOLDING_MODE_ALTERNATIVE
            shift
            ;;
        *)
            echo "ERROR: Unknown option $1"
            usage
            exit 1
            ;;
    esac
done

# check if both stage skip flags have been set
if [ "$SKIP_STAGE_1" = true ] && [ "$SKIP_STAGE_2" = true ]; then
    echo "Both --skip-stage-1 and --skip-stage-2 provided. Nothing to do. Exiting."
    exit 0
fi

if [ "$TEST_MODE" = true ]; then
    echo "
########################################
##### RUNNING FISHNET ON TEST DATA #####
########################################

    "

    ### test parameters ###
    trait="maleWC"
    traitRR="${trait}RR"
    num_permutations="10"
    pvalFileDir=./test/${trait}/
    pvalFileName=${pvalFileDir}/0-${trait}.csv
    pvalFileNameRR="./test/${traitRR}/${traitRR}.csv"
    module_algo="ker_based"
    moduleFileDir=./test/${module_algo}/
    numTests=$(( $(wc -l < "$pvalFileName") - 1 ))
    geneColName="Genes"
    pvalColName="p_vals"
    bonferroni_alpha="0.05"
    output_dir=./results/

    ### list of containers ###
    # contains all python dependencies for fishnet
    #   TODO: create single container with all python dependencies (include statsmodels)
    #   TODO: add to biocontainers 
    container_python="jungwooseok/dc_rp_genes:1.0"
    container_R="jungwooseok/r-webgestaltr:1.0"

    if [ "$SKIP_STAGE_1" = true ]; then
        echo "Skipping STAGE 1"
    else
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
        docker run --rm -v $(pwd):$(pwd) -w $(pwd) -u $(id -u):$(id -g) $container_python  /bin/bash -c \
            "python3 ./scripts/phase1/generate_uniform_pvals.py \
            --genes_filepath $pvalFileName"
        echo "done"


        # (3) nextflow random permutation run
        echo "# STEP 1.3: executing Nextflow MEA pipeline on random permutations"
        echo "executing Nextflow MEA pipeline on random permutations"
        nextflow run ./scripts/phase1/nextflow/main.nf \
            --trait ${traitRR} \
            --moduleFileDir $moduleFileDir \
            --numTests $numTests \
            --pipeline ${traitRR}  \
            --pvalFileName $pvalFileNameRR \
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
        summaries_path_original="${output_dir}/masterSummaries/summaries/"
        docker run --rm -v $(pwd):$(pwd) -w $(pwd) -u $(id -u):$(id -g) $container_python /bin/bash -c \
            "python3 ./scripts/phase1/compile_results.py \
                --dirPath $summaries_path_original \
                --identifier $trait \
                --output $output_dir"

        # (4.2) random permutation
        echo "# STEP 4.2: permutation results"
        summaries_path_permutation="${output_dir}/masterSummaries_RP/summaries/"
        docker run --rm -v $(pwd):$(pwd) -w $(pwd) -u $(id -u):$(id -g) $container_python /bin/bash -c \
            "python3 ./scripts/phase1/compile_results.py \
                --dirPath $summaries_path_permutation \
                --identifier ${traitRR} \
                --output $output_dir"
        echo "done"

        # (5) organize phase 1 results
        echo "# STEP 5: restructuring phase 1 results for input to phase 2"
        # TODO: consider defining shell functions in separate file
        # function to organize phase 1 nextflow pipeline results
        organize_nextflow_results() {
            trait=$1
            results_dir=$2
            traitRR="${traitRR}"

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
    fi

    if [ "$SKIP_STAGE_2" = true ]; then
        echo "Skipping STAGE 2"
    else
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

        ( head -n 1 "${output_dir}/${traitRR}/master_summary.csv"; grep 'True' "${output_dir}/${traitRR}/master_summary.csv" ) > "${output_dir}/${traitRR}/master_summary_filtered.csv"
        cut -d ',' -f 1-8 "${output_dir}/${traitRR}/master_summary_filtered.csv" > "${output_dir}/${traitRR}/master_summary_filtered_parsed.csv"
        
        if [ "$THRESHOLDING_MODE" = "$THRESHOLDING_MODE_DEFAULT" ]; then
            ##########################
            ## DEFAULT THRESHOLDING ##
            ##########################
            echo "
##########################
## DEFAULT THRESHOLDING ##
##########################
            "
            ## (1) generate statistics for original run
            echo "# STEP 2.1: generating statistics for original run"
            for network in `ls ${moduleFileDir}/`;
            do
                echo "Network: $network"
                network="${network%.*}" # remove extension
                docker run --rm -v $(pwd):$(pwd) -w $(pwd) -u $(id -u):$(id -g) $container_python /bin/bash -c \
                    "python3 ./scripts/phase2/dc_generate_or_statistics.py \
                        --gene_set_path $pvalFileDir \
                        --master_summary_path ${output_dir}/${trait}/master_summary_filtered_parsed.csv \
                        --trait $trait  \
                        --module_path ${moduleFileDir}/${network}.txt \
                        --go_path ${output_dir}/${trait}/GO_summaries/${trait}/ \
                        --study $trait \
                        --output_path ${output_dir}/${trait}/results/raw/ \
                        --network $network"
            done
            echo "done"

            # (2) generate statistics for permutation run
            # TODO: rewrite in nextflow(?) for parallelism
            # (2.1) TODO: prepare file with threshold:network pairs

            # (2.2) generate statistics for permutation run
            echo "# STEP 2.2: generating statistics for permutation runs"
            genes_rpscores_filedir="./results/RPscores/${traitRR}/"
            threshold_network_pairs=$( readlink -f ./test/slurm_thresholds_maleWC.txt )
            while IFS=$'\t' read -r threshold network; do
                echo "Threshold $threshold, Network: $network"
                docker run --rm -v $(pwd):$(pwd) -w $(pwd) -u $(id -u):$(id -g) $container_python /bin/bash -c \
                    "python3 ./scripts/phase2/dc_generate_rp_statistics.py \
                        --gene_set_path $genes_rpscores_filedir \
                        --master_summary_path ${output_dir}/${traitRR}/master_summary_filtered_parsed.csv \
                        --trait ${traitRR} \
                        --module_path ${moduleFileDir}/${network}.txt \
                        --go_path ${output_dir}/${traitRR}/GO_summaries/${traitRR}/ \
                        --output_path ${output_dir}/${traitRR}/results/raw/ \
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
                docker run --rm -v $(pwd):$(pwd) -w $(pwd) -u $(id -u):$(id -g) $container_python /bin/bash -c \
                    "python3 ./scripts/phase2/dc_summary_statistics_rp.py \
                        --trait $trait \
                        --input_path $output_dir \
                        --or_id $trait \
                        --rr_id ${traitRR} \
                        --input_file_rr_id ${traitRR} \
                        --network $network \
                        --output_path ${output_dir}/${trait}/summary/"
            done

            # (4) identify MEA passing genes
            echo "# STEP 2.4: identify MEA-passing genes"
            FDR_threshold=0.05 # TODO: add as input argument
            percentile_threshold=0.99 # TODO: add as input argument
            for network in `ls ${moduleFileDir}/`;
            do
                echo "Network: $network"
                network="${network%.*}" # remove extension
                docker run --rm -v $(pwd):$(pwd) -w $(pwd) -u $(id -u):$(id -g) $container_python /bin/bash -c \
                    "python3 ./scripts/phase2/dc_identify_mea_passing_genes.py \
                        --trait $trait \
                        --geneset_input $pvalFileDir\
                        --FDR_threshold $FDR_threshold \
                        --percentile_threshold $percentile_threshold \
                        --network $network \
                        --input_path ${output_dir}/${trait}"
            done
        else
            ##############################
            ## ALTERNATIVE THRESHOLDING ##
            ##############################
            echo "
##############################
## ALTERNATIVE THRESHOLDING ##
##############################
            "
            # (1) generate background gene sets for GO analysis
            echo "# STEP 1: generating background gene sets for GO analysis"
            docker run --rm -v $(pwd):$(pwd) -w $(pwd) -u $(id -u):$(id -g) $container_python /bin/bash -c \
                "python3 ./scripts/phase2/dc_fishnet_background_genes.py \
                    --genes_filepath $pvalFileNameRR \
                    --module_filepath $moduleFileDir \
                    --output_filepath ${output_dir}/${trait}/"

            # (1.1) copy background genes to permutation directory
            cp -r ${output_dir}/${trait}/background_genes/ ${output_dir}/${traitRR}/

            # (2) Extract and save module genes as individual files for modules that satisfy Bonferroni 0.25
            echo "# STEP 2: extract and save modules that satisfy Bonferroni threshold"
            # (2.1) original run
            echo " - original run"
            docker run --rm -v $(pwd):$(pwd) -w $(pwd) -u $(id -u):$(id -g) $container_python /bin/bash -c \
                "python3 ./scripts/phase2/dc_fishnet_module_genes.py \
                    --genes_filepath $pvalFileNameRR \
                    --module_filepath $moduleFileDir \
                    --master_summary_path ${output_dir}/${trait}/ \
                    --study $trait"
            # (2.1) permutation run
            echo " - permutation run"
            docker run --rm -v $(pwd):$(pwd) -w $(pwd) -u $(id -u):$(id -g) $container_python /bin/bash -c \
                "python3 ./scripts/phase2/dc_fishnet_module_genes.py \
                    --genes_filepath $pvalFileNameRR \
                    --module_filepath $moduleFileDir \
                    --master_summary_path ${output_dir}/${traitRR}/ \
                    --study $traitRR"

            # (3) Run GO analysis
            echo "# STEP 3: running GO analysis"
            # (3.1) original run
            echo " - original run"
            set +e
            for network in `ls ${moduleFileDir}/`;
            do
                echo "Network: $network"
                network="${network%.*}" # remove extension
                docker run --rm -v $(pwd):$(pwd) -w $(pwd) -u $(id -u):$(id -g) $container_R /bin/bash -c \
                    "Rscript ./scripts/phase2/dc_ORA_cmd.R \
                        --sigModuleDir ${output_dir}/${trait}/enriched_modules/${trait}-${network}/ \
                        --backGroundGenesFile ${output_dir}/${trait}/background_genes/${module_algo}-${network}.txt \
                        --summaryRoot ${output_dir}/${trait}/GO_summaries_alternate/ \
                        --reportRoot ${output_dir}/${trait}/report_sumamries_alternate/"

            done
            set -e
            # (3.2) permutation run
            echo " - permutation run"
            set +e
            for network in `ls ${moduleFileDir}/`;
            do
                echo "Network: $network"
                network="${network%.*}" # remove extension
                docker run --rm -v $(pwd):$(pwd) -w $(pwd) -u $(id -u):$(id -g) $container_R /bin/bash -c \
                    "Rscript ./scripts/phase2/dc_ORA_cmd.R \
                        --sigModuleDir ${output_dir}/${traitRR}/enriched_modules/${traitRR}-${network}/ \
                        --backGroundGenesFile ${output_dir}/${traitRR}/background_genes/${module_algo}-${network}.txt \
                        --summaryRoot ${output_dir}/${traitRR}/GO_summaries_alternate/ \
                        --reportRoot ${output_dir}/${traitRR}/report_sumamries_alternate/"

            done
            set -e

            # (4) Generate statistics for original gene-level p-values
            echo "# STEP 4: generate statistics for original gene-level p-values"
            for network in `ls ${moduleFileDir}/`;
            do
                echo "Network: $network"
                network="${network%.*}" # remove extension
                docker run --rm -v $(pwd):$(pwd) -w $(pwd) -u $(id -u):$(id -g) $container_python /bin/bash -c \
                    "python3 ./scripts/phase2/dc_generate_or_statistics_alternate.py \
                        --gene_set_path ${pvalFileDir} \
                        --master_summary_path ${output_dir}/${trait}/master_summary_alternate.csv \
                        --trait $trait \
                        --module_path ${output_dir}/${trait}/enriched_modules/${trait}-${network}/ \
                        --go_path ${output_dir}/${trait}/GO_summaries_alternate/ \
                        --study $trait \
                        --output_path ${output_dir}/${trait}/results/raw_alternate/ \
                        --network $network"
            done

            # (5) Generate statistics for random permutation runs
            echo "# STEP 5: generate statistics for random permutations"
            # (5.1) TODO: prepare file with bonferroni p-value:network pairs
            threshold_network_pairs_alternate=$( readlink -f ./test/slurm_thresholds_maleWC_alternate.txt )
            while IFS=$'\t' read -r threshold network; do
                echo "Threshold $threshold, Network: $network"
                docker run --rm -v $(pwd):$(pwd) -w $(pwd) -u $(id -u):$(id -g) $container_python /bin/bash -c \
                    "python3 ./scripts/phase2/dc_generate_rp_statistics_alternate.py \
                        --gene_set_path ${output_dir}/RPscores/${traitRR}/ \
                        --master_summary_path ${output_dir}/${traitRR}/master_summary_alternate.csv \
                        --trait ${traitRR} \
                        --module_path ${output_dir}/${traitRR}/enriched_modules/${traitRR}-${network}/ \
                        --go_path ${output_dir}/${traitRR}/GO_summaries_alternate/ \
                        --output_path ${output_dir}/${traitRR}/results/raw_alternate/ \
                        --network $network \
                        --threshold $threshold"
            done < $threshold_network_pairs_alternate
            echo "done"
        fi
    fi
else
    echo "FISHENT CURRENTLY ONLY SUPPORTS the --test FLAG"
fi
echo "### FISHNET COMPLETE ###"
