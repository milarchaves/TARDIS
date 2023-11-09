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
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from contrabass.core import CobraMetabolicModel

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

def find_essential_genes (map):
    '''Find essential genes in the metabolic map

    Parameters
    ----------
    model: CobraMetabolicModel
        Metabolic map in SBML format (.xml)

    Returns
    -------
    Essencial genes list

    Raises
    ------
    None
   '''
    
    model = CobraMetabolicModel(map)

    # update flux bounds with FVA
    model.fva(update_flux=True)

    # compute essencial genes
    model.find_essential_genes_1()

    # get essencial genes
    genes = model.essential_genes()
    
    if initial_args.verbosity > 0:
        print(clrs['g']+'Finding essential genes...'+clrs['n'])
        model.print_essential_genes()

    return genes
    

def find_chokepoint_reactions (map):
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
    
    model = CobraMetabolicModel(map)

    # update flux bounds with FVA
    model.fva(update_flux=True)

   # get chokepoints
    model.find_chokepoints(exclude_dead_reactions=True)
    chokepoints = model.chokepoints()

    chokepoint_list = []
    for i in range(len(chokepoints)):
        if chokepoints[i][0].id not in chokepoint_list:
            chokepoint_list.append(chokepoints[i][0].id)

    if initial_args.verbosity > 0:
        print(clrs['g']+'Finding chokepoint reactions...'+clrs['n'])
        model.print_chokepoints()

     # remove duplicates
    chokepoint_list = list(set(chokepoint_list))

    # get chokepoints
    return chokepoint_list

def find_essential_chokepoint_reactions (map):
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
    model = CobraMetabolicModel(map)

    # update flux bounds with FVA
    model.fva(update_flux=True)

    #find chokepoints
    model.find_chokepoints(exclude_dead_reactions=True)
    chokepoints = model.chokepoints()

    #find essential genes
    model.find_essential_genes_1()
    genes = list(model.essential_genes())
   
    # get essential_chokepoint_reactions
    essential_CP = []
    for i in range(len(chokepoints)):
        for j in range(len(genes)):
            if chokepoints[i][0].gene_reaction_rule == genes[j].id:
                essential_CP.append(genes[j].id)

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
  

    
