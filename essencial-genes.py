'''
++++++ Imports ++++++
'''
import pandas as pd
import cobra
from cobra.flux_analysis import ( single_gene_deletion, single_reaction_deletion, double_gene_deletion, double_reaction_deletion)
from cobra.flux_analysis import flux_variability_analysis

'''
++++++ File path ++++++
'''
model = '/home/camila/LMDM/P.Aeruginosa/CCBH4851/CCBH4851-metabolic-network/CCBH4851_metabolic_map_v2.xml'

#open sbml metabolic model
model = cobra.io.read_sbml_model(model)

#maximize flux through the objective reactions
model.optimize()

#finds the ranges of each metabolic flux at 90% optimality (literature review)
model.summary(fva=0.9)

#perform all single gene deletions on a model
deletion = single_gene_deletion(model)

#finds the genes wich deletion impacts growth rate
impact_genes = deletion[deletion['growth'] < 12.470972]

#finds the genes wich deletion abolish growth rate
essencial_genes = deletion[deletion['growth'] <= 0]

essencial_genes.to_csv('essen_genes.csv')

df = pd.read_csv('essen_genes.csv')
print(df['ids'])
