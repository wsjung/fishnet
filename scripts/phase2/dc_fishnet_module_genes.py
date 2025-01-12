import numpy as np
import pandas as pd
import ast
import os
import math
import scipy.stats as stats


# python dc_fishnet_module_genes.py 
#   --genes_filepath /scratch/mblab/acharyas/fishnet/pipeline/data/pvals/llfsTWASRR/llfsTWASRR.csv 
#   --module_filepath /scratch/mblab/acharyas/fishnet/pipeline/data/modules/ker_based/ 
#   --master_summary_path /scratch/mblab/acharyas/fishnet/pipeline/results/twasLLFSORKB/ 
#   --study twasLLFSORKB


def dc_fishnet_module_genes(genes_filepath, module_filepath, master_summary_path, study ):
    module_network_pair = extract_modules(master_summary_path, study)
    extract_module_genes(study, genes_filepath, module_filepath, module_network_pair, master_summary_path)

def extract_modules(master_summary_path, study):    
    final_module_network_pair = pd.DataFrame(columns = ["moduleIndex", "network", "trait", "Study"])
    master_summary_df = pd.read_csv(os.path.join(master_summary_path,"/master_summary.csv"))
    master_summary_df = master_summary_df[["study", "trait", "network", "moduleIndex", "isModuleSig", "modulePval", "moduleBonPval", "size"]]
    master_summary_df = master_summary_df[master_summary_df["moduleBonPval"] <= 0.25]
    master_summary_df.to_csv(os.path.join(master_summary_path,"master_summary_alternate.csv"), index = False)

    final_module_network_pair = pd.DataFrame(columns = ["moduleIndex", "network", "trait", "Study"])

    for networks in master_summary_df["network"].unique():
        for traits in master_summary_df["trait"].unique():
            temp_master_summary_df = master_summary_df[(master_summary_df["trait"] == traits) &
                                                        (master_summary_df["network"] == networks)]
            temp_master_summary_df = temp_master_summary_df[["moduleIndex", "network", "trait"]]
            temp_master_summary_df["Study"] = study
            final_module_network_pair = pd.concat([final_module_network_pair,temp_master_summary_df], ignore_index = True )
    final_module_network_pair.to_csv(os.path.join(master_summary_path,f"{study}.txt"), index = False)
    final_module_network_pair = final_module_network_pair[["moduleIndex", "network"]]
    filtered_final_module_network_pair = final_module_network_pair.drop_duplicates(subset=['moduleIndex', 'network'], ignore_index = True)
    return filtered_final_module_network_pair
 
def save_module_genes(modules_df, study, genes_filepath, master_summary_path):
    for index, rows in modules_df.iterrows():
        genes = rows["Genes"]
        module_id = rows["ID"]
        network = rows["Network_type"]
        module_genes = genes.split("\t")
        
        #take intersection of module genes and summary stats genes
        gene_summary_df = pd.read_csv(genes_filepath)
        genes_summary = gene_summary_df["Genes"]
        final_genes = set(genes_summary).intersection(set(module_genes))
        output_dir = os.path.join(master_summary_path,"enriched_modules",f"{study}-{network}")
        os.makedirs(output_dir, exist_ok=True)

        with open(os.path.join(output_dir,f"sig_{study}-{network}-{module_id}.txt"), 'w') as file:
            for item in genes.split("\t"):
                file.write(f"{item}\n") 

def extract_module_genes(study, genes_filepath, module_filepath, module_network_pair, master_summary_path):
    Index = []
    Network = []
    Genes = []

    for files in os.listdir(module_filepath):
        network_modules = module_network_pair[module_network_pair["network"] == files.split(".")[0]]["moduleIndex"]
        if(network_modules.shape[0] == 0):
            continue
        with open(os.path.join(module_filepath,files), 'r') as file:
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
        modules_df = pd.DataFrame({
            'ID': Index,
            'Network': Network,
            'Genes': Genes
        })
        modules_df["ID"] = modules_df["ID"].astype(int)
        modules_df = modules_df[modules_df["ID"].isin(network_modules)]
        modules_df.loc[:, "Network_type"] = files.split(".")[0]        
        save_module_genes(modules_df, study, genes_filepath, master_summary_path)


if __name__ == "__main__":
    from argparse import ArgumentParser   
    parser = ArgumentParser()
    parser.add_argument('--genes_filepath', '-genes', help='genes')
    parser.add_argument('--module_filepath', '-module_filepath')
    parser.add_argument('--master_summary_path', '-master_summary_path')
    parser.add_argument('--study', '-study')

    args = parser.parse_args()
    dc_fishnet_module_genes(genes_filepath = args.genes_filepath, module_filepath = args.module_filepath, master_summary_path =  args.master_summary_path, study = args.study)


