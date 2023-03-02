'''
++++++ Imports ++++++
'''
import pandas as pd
import cobra
from findCPcore import CobraMetabolicModel
from cobra.flux_analysis import ( single_gene_deletion, single_reaction_deletion, double_gene_deletion, double_reaction_deletion)
from cobra.flux_analysis import flux_variability_analysis

'''
++++++ File path ++++++
'''
model = '/home/camila/LMDM/Mestrado/TARDIS/test/d_radiodurans.xml'

#open sbml metabolic model
model = cobra.io.read_sbml_model(model)

model = CobraMetabolicModel('/home/camila/LMDM/Mestrado/TARDIS/test/d_radiodurans.xml')

#maximize flux through the objective reactions
model.optimize()

# update flux bounds with FVA
model.fva(update_flux=True)

# compute chokepoints
model.find_chokepoints()


# get chokepoints
model.find_chokepoints(exclude_dead_reactions=True)
chokepoints = model.chokepoints()
chokepoint_list = []
for i in range(len(chokepoints)):
    reaction = chokepoint[i].__str__().split(' ')
    chokepoint_list.append(reaction[1])

chokepoint_list

dir(chokepoints)
model.find_essential_genes_reactions()

model.essential_genes_reactions()

model.find_essential_genes_1()

model.essential_genes()


#Generates chokepoint summary spreadsheet
from findCPcore import Facade
from findCPcore import FacadeUtils

def callback_print_ignore(message, arg1, arg2):
    pass

def callback_print_logger(message, arg1, arg2):
    logger = arg1
    logger.info(message)

def callback_print(message, arg1, arg2):
    print(arg1 + message)

facadeUtils = FacadeUtils()

dir(FacadeUtils())
spreadsheet = facadeUtils.run_summary_model('/home/camila/LMDM/Mestrado/TARDIS/test/d_radiodurans.xml', callback_print, "LOG:", None, fraction=0.95)
facadeUtils.save_spreadsheet("output.xls", spreadsheet)


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
