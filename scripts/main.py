import subprocess
import argparse
import pandas as pd
import sys

from dc_rp_genes import dc_rp_genes

def nextflow_cleanup():
    """This function cleans up the nextflow work directory and logs"""

    print("cleaning up nextflow logs and work/")

    # clean up work/
    subprocess.run(["nextflow", "clean", "-f"])
    # delete nextflow log files
    subprocess.run(["rm", "-rf", ".nextflow/", ".nextflow.log*", "work/"])

    print("done")


def run_nexflow_mea_original(trait, input_file, module_dir_path, debug=False):
    """This function runs the nextflow module and GO enrichment analysis pipeline

    Parameters
    ----------
    trait : string
            Trait name
    input_file : string
            Path to the input file
    module_dir_path : string
            Path to the directory containing module files
    """

    pval_file = pd.read_csv(input_file, sep=",")
    num_tests = pval_file.shape[0]

    params = [
        "--trait", trait,
        "--moduleFileDir", module_dir_path,
        "--geneColName", "Genes",
        "--pvalColName", "p_vals",
        "--numTests", str(num_tests),
        "--pipeline", trait,
        "--pvalFileName", input_file
    ]

    if debug:
        print(["nextflow", "run", "/app/scripts/scripts_nf/mea_slurm_or_kb.nf"] + params)

    subprocess.run(["nextflow", "run", "/app/scripts/scripts_nf/mea_slurm_or_kb.nf"] + params, check=True)

    if not debug:
        nextflow_cleanup()

def generate_uniformly_distributed_pvals(input_file):
    """This function generates uniformly distributed p-values

    Storse the generated p-values in a new directory with
    input file directory name appended with RR

    Parameters
    ----------
    input_file : string
            Path to the input file

    Returns
    -------
    new_file_path : string
            Path to the generated p-values file
    """
    new_file_path = dc_rp_genes(input_file)
    return new_file_path

def run_nextflow_mea_random_permutation(trait, input_file, module_dir_path,
                                        num_permutations, debug=False):
    """This function runs the nextflow module and GO enrichment analysis pipeline

    Parameters
    ----------
    trait : string
            Trait name
    input_file : string
            Path to the input file with generated p-values
    module_dir_path : string
            Path to the directory containing module files
    num_permutations : int
            Number of random permutation runs
    debug : bool
            Set to True to turn on debugging (default: False)
    """
    pval_file = pd.read_csv(input_file, sep=",")
    num_tests = pval_file.shape[0]

    params = [
        "--trait", trait,
        "--moduleFileDir", module_dir_path,
        "--geneColName", "Genes",
        "--pvalColName", "p_vals",
        "--numTests", str(num_tests),
        "--pipeline", trait,
        "--pvalFileName", input_file,
        "--numRP", str(num_permutations)
    ]

    if debug:
        print(["nextflow", "run", 
               "/app/scripts/scripts_nf/mea_slurm_rp.nf"] + params)

    subprocess.run(["nextflow", "run",
                    "/app/scripts/scripts_nf/mea_slurm_rp.nf"] + params,
                   check=True)

    if not debug:
        nextflow_cleanup()


def run_test(debug):
    trait = "maleWC"
    trait_file_path = "/app/test/0-maleWC.csv"
    module_dir_path = "/app/test/ker_based/"
    num_permutations = 3
    #run_nexflow_mea_original(trait, trait_file_path, module_dir_path, debug)
    trait_RR_file_path = generate_uniformly_distributed_pvals(trait_file_path)
    run_nextflow_mea_random_permutation(trait, trait_RR_file_path,
                                        module_dir_path, num_permutations, debug)


def main():
    parser = argparse.ArgumentParser(description="help")
    parser.add_argument("--trait", required=False, help="trait")
    parser.add_argument("--trait_file_path", required=False, help="path to trait file")
    parser.add_argument("--module_directory_path", required=False, help="path to module directory")
    parser.add_argument("--test", action='store_true', help="run with minimal test data")
    parser.add_argument("--debug", action='store_true', help='Use to turn on debugging (will not delete the nextflow work/ directories and nextflow logs')
    args = parser.parse_args()

    if args.test or not len(sys.argv) > 1:
        print("initiating test run")
        run_test(args.debug)
        print("test run complete")
        raise SystemExit

    print("running nextflow")
    run_nexflow_mea_original(args.trait_file_path, args.trait, args.module_directory_path, args.debug)
    print("complete")

    print("generating pvals")
    trait_RR_file_path = generate_uniformly_distributed_pvals(args.trait_file_path)
    print("complete")


if __name__ == '__main__':
    main()
