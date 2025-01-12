import numpy as np
import pandas as pd
import ast
import os
import math
import scipy.stats as stats


# python dc_fishnet_background_genes.py 
#   --genes_filepath /scratch/mblab/acharyas/fishnet/pipeline/data/pvals/llfsTWASRR/llfsTWASRR.csv 
#   --module_filepath /scratch/mblab/acharyas/fishnet/pipeline/data/modules/ker_based/ 
#   --output_filepath /scratch/mblab/acharyas/fishnet/pipeline/results/twasLLFSORKB/

def dc_fishnet_background_genes(genes_filepath, module_filepath, output_filepath ):
    networks = []
    for files in os.listdir(module_filepath):
        genes = background_set(os.path.join(module_filepath,files), genes_filepath)
        moduleAlgo = os.path.normpath(module_filepath)
        moduleAlgo = os.path.basename(moduleAlgo)
        output_dir = os.path.join(output_filepath,"background_genes")
        os.makedirs(output_dir, exist_ok = True)
        with open(os.path.join(output_dir,f"{moduleAlgo}-{files}"), 'w') as file:
            # Write each gene to a new line
            for gene in genes:
                file.write(f"{gene}\n")

def background_set(network_path, genes_filepath):
    Index = []
    Network = []
    Genes = []
    #coex
    with open(network_path, 'r') as file:
            for line in file:
                # Split the line into components
                parts = line.strip().split()
                # The first element is the ID
                Index.append(parts[0])

                # The second element is the Attribute
                Network.append(parts[1])

                # The rest are genes
                Genes.append('\t'.join(parts[2:]))

            # Create a DataFrame from the lists
    
    modules_df= pd.DataFrame({
        'ID': Index,
        'Network': Network,
        'Genes': Genes
    })
    genes_list = modules_df['Genes'].str.split('\t').sum()
    genes_summary_df = pd.read_csv(genes_filepath)
    genes_summary_list = genes_summary_df["Genes"]

    intersecting_genes = set(genes_summary_list).intersection(genes_list)
    print(len(intersecting_genes))
    return intersecting_genes    

def get_key_from_value(dictionary, value):
    for key, val in dictionary.items():
        if val == value:
            return key
    return None  # Return None if the value is not found

 
if __name__ == "__main__":
    from argparse import ArgumentParser   
    parser = ArgumentParser()
    parser.add_argument('--genes_filepath', '-genes', help='genes')
    parser.add_argument('--module_filepath', '-module_filepath')
    parser.add_argument('--output_filepath', '-output_filepath')
    args = parser.parse_args()
    dc_fishnet_background_genes(genes_filepath = args.genes_filepath, module_filepath = args.module_filepath, output_filepath = args.output_filepath)


