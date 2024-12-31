import subprocess
import argparse
import sys

def run_nextflow(trait, input_file, module_dir_path):

    params = [
        "--trait", trait,
        "--moduleFileDir", module_dir_path,
        "--geneColName", "Genes",
        "--pvalColName", "p_vals",
        "--numTests",
        "--pipeline","maleWC",
        "--pvalFileName", input_file
    ]

    subprocess.run(["nextflow", "run", "/app/scripts/mea_slurm_or_kb.nf"] + params, check=True)

def run_test():
    trait = "maleWC"
    trait_file_path = "/app/test/maleWC.csv"
    module_dir_path = "/app/test/ker_based/"
    run_nextflow(trait, trait_file_path, module_dir_path)

def main():
    parser = argparse.ArgumentParser(description="help")
    parser.add_argument("--trait", required=False, help="trait")
    parser.add_argument("--trait_file_path", required=False, help="path to trait file")
    parser.add_argument("--module_directory_path", required=False, help="path to module directory")
    parser.add_argument("--test", action='store_true', help="run with minimal test data")
    args = parser.parse_args()

    if args.test or not len(sys.argv) > 1:
        print("initiating test run")
        run_test()
        print("test run complete")
        raise SystemExit

    run_nextflow(args.trait_file_path, args.trait, args.module_directory_path)


if __name__ == '__main__':
    main()
