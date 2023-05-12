#Imports
###############################################################################
import pandas as pd

# Paths
deg = "/home/camila/LMDM/Mestrado/TARDIS/deg_all_bacteria/deg_annotation_p.csv"
bacillus_tardis = "/home/camila/LMDM/Mestrado/TARDIS/essential_genes_validation/bacillus_tardis.txt"
gene_names = "/home/camila/LMDM/Mestrado/TARDIS/essential_genes_validation/bacillus_tardis_gene_name.txt"



#open deg file as a pandas dataframe
df = pd.read_csv(deg, sep = ",", index_col=False, names=["Gene", "Function", "Description", "Bacteria", "Uniprot"])

#Get only Bacillus subtilis 168 data
df_bacillus = df[df['Bacteria']== "Bacillus subtilis 168"]


#open bacillus_tardis file remove \n and split by "_"
essential_genes = []
with open(bacillus_tardis, 'r') as f:
    for line in f:
        essential_genes.append(line.split("_")[1])

essential_genes = [x.replace('\n', '') for x in essential_genes]

#compare the two lists
tardis_deg = []
for gene in essential_genes:
    if gene in df_bacillus["Uniprot"].values:
        tardis_deg.append(gene)
        print(gene)

#save the results in a file
with open('tardis_deg.txt', 'w') as f:
    for item in tardis_deg:
        f.write("%s\n" % item)

#save genes from bacillus_tardis that are not in deg file
tardis_not_deg = []
for gene in essential_genes:
    if gene not in df_bacillus["Uniprot"].values:
        tardis_not_deg.append(gene)
print(len(tardis_not_deg))

#save the results in a file
with open('tardis_not_deg.txt', 'w') as f:
    for item in tardis_not_deg:
        f.write("%s\n" % item)

#Count the number of genes in each list
print(len(essential_genes))
print(len(tardis_deg))
print(len(df_bacillus["Uniprot"].values))

#open gene_names file remove \n and split by " "

tardis_gene_names = []
with open(gene_names, 'r') as f:
    for line in f:
        tardis_gene_names.append(line.split(" ")[0])

print(tardis_gene_names)
print(df_bacillus["Gene"])

#compare the two lists
tardis_deg_gene_names = []
for gene in tardis_gene_names:
    if gene in df_bacillus["Gene"].values:
        tardis_deg_gene_names.append(gene)
        print(gene)
print(len(tardis_deg_gene_names))