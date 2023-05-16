#!/usr/lib/python3

# Description
###############################################################################
'''
This file contain the metabolic moddeling routine to find new targets.

1. Create the metabolic map using CarveMe (D. Machado et al, 2018. https://doi.org/10.1093/nar/gky537)
2. Find essential genes
3. Find chokepoint reactions (Oarga et al, 2020. https://doi.org/10.1007/978-3-030-60327-4_6)
4. Check if exists homology in human genome

A protein is considered a potential target if it fits all the criteria above.  

'''

# Imports
###############################################################################
from TARDIS.Initialise import *
import subprocess
import cobra
from contrabass.core import CobraMetabolicModel
from cobra.flux_analysis import single_gene_deletion

# License
###############################################################################
'''

TARDIS: TARgets DIScoverer

Authors: Chaves, C; Torres, P.H.M.

[Federal University of Rio de Janeiro]

Contact info:
E-mail address: chaves.camila13@gmail.com
Github: https://github.com/milarchaves
This project is licensed under Creative Commons license (CC-BY-4.0)

'''

# functions
###############################################################################

def create_map (input):
    '''Open input file, parse protein name and sequence to output dataframe and create the metabolic map

    Parameters
    ----------
    Imput file in FASTA format (.faa)

    Returns
    -------
    Metabolic map in SBML format (.xml)

    Raises
    ------
    No input file or input file in wrong format error
    '''
    if input is None:
        print(clrs['r']+'ERROR: '+clrs['n']+'No input file was provided. Please use the'+clrs['y']+' -f '+clrs['n']+'flag to provide the proteome file.')
    if input.endswith('.faa') or input.endswith('.fasta'):
        if initial_args.verbosity > 0:
            print(clrs['g']+'Input file provided. Creating metabolic map...'+clrs['n'])
        subprocess.run(["carve", input])
        if input.endswith('.faa'):
            model = input.replace('.faa', '.xml')
        else:
            model = input.replace('.fasta', '.xml')
    else:
        print(clrs['r']+'ERROR: '+clrs['n']+'The input file must be in FASTA format (.faa).')
    return model

def refine_model (input):
    '''Refine the metabolic map

    Parameters
    ----------
    Metabolic map in SBML format (.xml)

    Returns
    -------
    Refined metabolic map in SBML format (.xml)
    '''

    #open sbml metabolic model
    model = cobra.io.read_sbml_model(input)

    #Open the model as a CONSTRABASS object to use the methods below
    CBmodel = CobraMetabolicModel(model)

    # update flux bounds with FVA
    CBmodel.fva(update_flux=True)

    #remove Dead-end metabolites
    CBmodel.remove_dem()

    #save the refined model
    cobra.write_sbml_model(CBmodel, (input.replace('.xml', '_refined.xml')))
    refined_model = input.replace('.xml', '_refined.xml')

    return refined_model




def find_essential_genes (model):
    '''Find essential genes in the metabolic map

    Parameters
    ----------
    Metabolic map in SBML format (.xml)

    Returns
    -------
    Essencial genes list

    Raises
    ------
    None
   '''

    #open sbml metabolic model
    model = cobra.io.read_sbml_model(model)

    #open sbml metabolic model
    model = cobra.io.read_sbml_model('/home/camila/LMDM/Mestrado/TARDIS/TARDIS/Vibrio_cholerae_serotype_O1.xml')

    #Open the model as a CONSTRABASS object to use the methods below
    CBmodel = CobraMetabolicModel('cami/home/la/LMDM/Mestrado/TARDIS/TARDIS/Vibrio_cholerae_serotype_O1.xml')

    # update flux bounds with FVA
    CBmodel.fva(update_flux=True)

    #remove Dead-end metabolites
    CBmodel.remove_dem()

    dir(CBmodel)

    #maximize flux through the objective reactions
    CBmodel.optimize()

    #Performe FVA
    model.summary(fva=0.9)

    #perform all single gene deletions on a model
    solution = CBmodel.optimize()
    solution.objective_value

    #perform all single gene deletions on a model
    deletion = single_gene_deletion(model)
    deletion = deletion.reset_index()

    #Growth trashold
    trashold = 0.25*float(solution.objective_value)
    trashold

    #finds genes which, when deleted, result in growth less than 5% of the wild-type growth
    essencial_genes = deletion[deletion['growth'] <= trashold]
    non_essencial_genes = deletion[deletion['growth'] > trashold]
    print(non_essencial_genes)


    #essencial_genes = deletion[deletion['growth'] <= (0.05*solution)]

    return essencial_genes['ids']
    

def find_chokepoint_reactions (model):
    '''Find chokepoint reactions in the metabolic map

    Parameters
    ----------
    Metabolic map in SBML format (.xml)

    Returns
    -------
    Chokepoint reactions list

    Raises
    ------
    None
   '''
         
    #Open model
    cobramodel = CobraMetabolicModel(model)

    # update flux bounds with FVA
    cobramodel.fva(update_flux=True)

   # get chokepoints
    cobramodel.find_chokepoints(exclude_dead_reactions=True)
    chokepoints = cobramodel.chokepoints()

    chokepoint_list = []
    for i in range(len(chokepoints)):
        if chokepoints[i][0].id not in chokepoint_list:
            chokepoint_list.append(chokepoints[i][0].id)

    if initial_args.verbosity > 0:
        print(clrs['g']+'Finding chokepoint reactions...'+clrs['n'])
        cobramodel.print_chokepoints()

     # remove duplicates
    chokepoint_list = list(set(chokepoint_list))

    # get chokepoints
    return chokepoint_list

def find_essential_chokepoint_reactions (model, essential_genes):
    '''Find essential genes reactions in the metabolic map

    Parameters
    ----------
    Metabolic map in SBML format (.xml)

    Returns
    -------
    essential genes reactions list

    Raises
    ------
    None
   '''
    
    #Open model
    model = CobraMetabolicModel(model)
    #find chokepoints
    model.find_chokepoints(exclude_dead_reactions=True)
    chokepoints = model.chokepoints()

    #find essential genes
    genes = list(essential_genes)
   
    # get essential_chokepoint_reactions
    essential_CP = []
    for i in range(len(chokepoints)):
        for j in range(len(genes)):
            if chokepoints[i][0].gene_reaction_rule == genes[j]:
                essential_CP.append(genes[j])

    # remove duplicates
    essential_CP = list(set(essential_CP))

    if initial_args.verbosity > 0:
        print(clrs['g']+'Essential chokepoint genes...'+clrs['n'])
        print(essential_CP)
   
    return essential_CP

def homology_search (model, output):
    '''Find homology in human genome

    Parameters
    ----------
    Metabolic map in SBML format (.xml)
    Output dataframe

    Returns
    -------
    Output dataframe with homology information

    Raises
    ------
    None
   '''
    if initial_args.verbosity > 0:
        print(clrs['g']+'Finding homology in human genome...'+clrs['n'])
    #find homology in human genome
  

    
