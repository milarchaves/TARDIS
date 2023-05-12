#!/usr/lib/python3

# Imports
###############################################################################

import os
import pandas as pd
from matplotlib_venn import venn2
import matplotlib.pyplot as plt

# Paths
base_dir = "/home/camila/LMDM/Mestrado/TARDIS/essential_genes_validation/"
deg = "/home/camila/LMDM/Mestrado/TARDIS/deg_all_bacteria/deg_annotation_p.csv"

# Functions
###############################################################################

def essential_gene_validation(file):
 
    #open deg file as a pandas dataframe
    df = pd.read_csv(deg, sep = ",", index_col=False, names=["Gene", "Function", "Description", "Bacteria", "Uniprot"])

    #Remove duplicates and NA from df[Bacteria]
    df['Bacteria'] = df['Bacteria'].fillna('NA')
    df = df.drop_duplicates(subset=['Bacteria'])
    print(df['Bacteria'])

    #open file remove \n and split by "_"
    tardis_eg = []
    with open(file, 'r') as f:
        for line in f:
            tardis_eg.append(line.split("_")[1])
    #parse file name replacing "_" by " " and removing ".txt"
    file_name = file.replace("_", " ").replace(".txt", "")

    if file_name in df['Bacteria'].values:
        spec_df = df[df['Bacteria']== file_name]
        #compare the two lists
        tardis_deg = []
        for gene in tardis_eg:
            if gene in spec_df["Uniprot"].values:
                tardis_deg.append(gene)
                print(gene)
        #save the results in a file
        with open(file.replace(".txt", "") + '_deg.txt', 'w') as f:
            for item in tardis_deg:
                f.write("%s\n" % item)
        #save genes from tardis_eg that are not in spec_df["Uniprot"]
        tardis_not_deg = []
        for gene in tardis_eg:
            if gene not in spec_df["Uniprot"].values:
                tardis_not_deg.append(gene)
        with open(file.replace(".txt", "") + '_not_deg.txt', 'w') as f:
            for item in tardis_not_deg:
                f.write("%s\n" % item)
        #save genes from spec_df["Uniprot"] that are not in tardis_eg
        deg_not_tardis = []
        for gene in spec_df["Uniprot"].values:
            if gene not in tardis_eg:
                deg_not_tardis.append(gene)
        with open(file.replace(".txt", "") + '_not_tardis.txt', 'w') as f:
            f.write("%s\n" % item)
        #print quantities
        print("Total essential genes in TARDIS: " + str(len(tardis_eg)))
        print("Total essential genes in DEG: " + str(len(spec_df["Uniprot"].values)))
        print("Essential genes in TARDIS and DEG: " + str(len(tardis_deg)))
        print("Essential genes in TARDIS and not in DEG: " + str(len(tardis_not_deg)))
        print("Essential genes in DEG and not in TARDIS: " + str(len(deg_not_tardis)))
        #plot venn diagram in pastel colors with transparency high quality and file name as title and save it in png
        plt.figure(figsize=(4,4))
        venn2(subsets = (len(tardis_not_deg), len(deg_not_tardis), len(tardis_deg)), set_labels = ('TARDIS', 'DEG'), set_colors=('skyblue', 'pink'), alpha = 0.5)
        plt.title(file_name)
        plt.savefig(file.replace(".txt", "") + '.png', dpi=300)
        plt.show()
    else:
        print("Bacteria not found in deg file")


# Go to basedir and run essential genes validation for each file
os.chdir(base_dir)
for file in os.listdir(base_dir):
    essential_gene_validation(file)
