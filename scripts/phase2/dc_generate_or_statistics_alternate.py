import numpy as np
import pandas as pd
import ast
import os
import math

#Usage in generate_or_statistics.sbatch

def generate_or_statistics(gene_set_path, master_summary_path, trait, module_path, go_path, study, output_path, network ):

    print("Started")
    final_df = pd.DataFrame(columns = ["threshold", "mea_passing_genes"])
    mea_passing_genes_df = pd.DataFrame(columns = ["threshold", "mea_passing_genes"])

    #change trait name to reflect the name of the input file
    trait = "0-" + trait


    #read gene sets df
    gene_set_df = pd.read_csv(os.path.join(gene_set_path,f"{trait}.csv"))
    gene_set_df.columns = ["Gene", "pval"]
    gene_set_df = gene_set_df.sort_values(by = ["pval"])

    #pick threshold to include top 5% genes.
    gene_set_threshold = int(gene_set_df.shape[0] * 0.05)

    temp_gene_set_df = gene_set_df.head(gene_set_threshold)
    gene_set_series = temp_gene_set_df["Gene"]


    #threshold for module enrichment analysis bon-corrected p-value%. 
    thresholds_range = [0.00005, 0.0001, 0.005, 0.01, 0.05, 0.1, 0.15, 0.20, 0.25]
    
    #read and process master summary file
    master_summary = pd.read_csv(master_summary_path)
    temp_master_summary = master_summary[master_summary["trait"] == trait]
    if (temp_master_summary.shape[0]) > 0:
        temp_master_summary = temp_master_summary[(temp_master_summary["network"] == network) &          
                                   (temp_master_summary["study"] == study)] 
        
        temp_master_summary = temp_master_summary.reset_index(drop=True)
    else:
        print("no enriched module for the trait: " + trait)
        return

    #iterate through multiple threshold of module enrichment significant p value
    for threshold in thresholds_range: 

        #get the number of MEA passing genes from the queried gene set
        mea_output = MEA_passing(module_path, temp_master_summary, go_path, gene_set_series, trait, study, network, threshold)
        temp_mea_passing_genes_count = mea_output[0]
        temp_mea_passing_genes = mea_output[1]        
        final_df.loc[len(final_df.index)] = [threshold, temp_mea_passing_genes_count]
        mea_passing_genes_df.loc[len(mea_passing_genes_df.index)] = [threshold, list(temp_mea_passing_genes) ]

    #write final df to the output directory
    final_df.to_csv(os.path.join(output_path,f"{network}_{trait}_{network}_or_summary.csv"), index = None)
    mea_passing_genes_df.to_csv(os.path.join(output_path,f"{network}_{trait}_{network}_or_fishnet_genes.csv"), index = None)

    print("Completed")

def get_GO_genes(go_df):
    go_df = go_df[["geneSet", "description", "size", "overlap", "FDR", "userId"]]
    go_df = go_df.sort_values(by = ["FDR"])
    
    #FDR threshold is set to 0.05
    go_df = go_df[go_df["FDR"] <= 0.05]

    #filter for genes annoted in enriched GO Term 
    all_go_genes = []
    for index, row in go_df.iterrows():
        temp_go_set = row["userId"].split(";")
        all_go_genes.extend(temp_go_set)
    all_go_genes_set = set(all_go_genes)
    return all_go_genes_set


def MEA_passing(module_path, temp_master_summary, go_path, gene_set, trait, study, network, threshold):
        mea_passing_genes = []
        temp_master_summary = temp_master_summary[temp_master_summary["moduleBonPval"] <= threshold]
        for index in temp_master_summary.index:
            #find the genes that lie in enriched modules
            module_index = temp_master_summary["moduleIndex"][index]
            module_index = int(module_index)
            module_df = pd.read_csv(os.path.join(module_path,f"sig_{module_path.split("/")[-2]}-{module_index}.txt"), header = None)
            module_df.columns = ["Genes"]
            module_genes = module_df["Genes"].tolist()
     
            #genes intersecting between the first criteria and the second criteria for MEA-passing genes for a given module
            module_genes_intersection_gene_set = set(gene_set).intersection(set(module_genes))

            #find the genes that also lie in enriched GO Term for the enriched module        
            go_file = f"sig_{module_path.split("/")[-2]}-{module_index}.csv"

            go_df = pd.read_csv(os.path.join(go_path,go_file))
            if (go_df.shape[0] > 0):
                all_go_genes_set = get_GO_genes(go_df)
                set_qualifying_genes = module_genes_intersection_gene_set.intersection(all_go_genes_set)
                mea_passing_genes.extend(list(set_qualifying_genes))

        #return the count of the number of MEA_passing genes
        set_mea_passing_genes = set(mea_passing_genes)
        return (len(set_mea_passing_genes), mea_passing_genes)

if __name__ == "__main__":
    from argparse import ArgumentParser   
    parser = ArgumentParser()
    parser.add_argument('--gene_set_path', '-gene_set_path', help='the path to the file that has genes and pvalues for a given trait')
    parser.add_argument('--master_summary_path', '-master_summary_path', help='the path to the master summary file')
    parser.add_argument('--trait', '-trait', help='trait')
    parser.add_argument('--module_path', '-module_path', help='path to the module genes')
    parser.add_argument('--go_path', '-go_path', help='path to enriched go terms')
    parser.add_argument('--study', '-study', help='study')
    parser.add_argument('--output_path', '-output_path', help = "directory to store the parsed MMAP output")
    parser.add_argument('--network', '-network', help = "network type")
    args = parser.parse_args()
    generate_or_statistics(gene_set_path = args.gene_set_path, master_summary_path = args.master_summary_path, trait = args.trait, module_path = args.module_path, go_path = args.go_path, study = args.study, output_path = args.output_path, network = args.network)


