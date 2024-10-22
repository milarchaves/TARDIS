#!/usr/bin/env python3

# Description
###############################################################################
'''
Routines to find new targets in the metabolic map.

1. Find essential genes
2. Find chokepoint reactions (Oarga et al, 2020. https://doi.org/10.1007/978-3-030-60327-4_6)
3. Check if exists homology in human genome

A protein is considered a potential target if it fits all the criteria above.

They are imported as:

import TARDIS.FindTargets as FindTargets
'''

# Imports
###############################################################################
from TARDIS.Initialise import *

from contrabass.core import CobraMetabolicModel
from cobra import Gene

# License
###############################################################################
'''
TARDIS: TARgets DIScoverer

Authors: Chaves, C; Rossi, A.D; Torres, P.H.M.

[Federal University of Rio de Janeiro]

Contact info:
E-mail address: chaves.camila13@gmail.com
Github: https://github.com/milarchaves
This project is licensed under Creative Commons license (CC-BY-4.0)
'''

# functions
###############################################################################

def find_essential_genes(map: str) -> set[Gene]:
    '''Find essential genes in the metabolic map

    Parameters
    ----------
    map (str)
        Path to the metabolic map in SBML format (.xml)

    Returns
    -------
    set[cobra.Gene]
        Essencial genes list
    '''
    
    model = CobraMetabolicModel(map)

    # update flux bounds with FVA
    model.fva(update_flux = True)

    # compute essencial genes
    model.find_essential_genes_1()

    # get essencial genes
    genes = model.essential_genes()
    
    if initial_args.verbosity > 0:
        print(f"{clrs['g']}Finding essential genes...{clrs['n']}")
        model.print_essential_genes()

    return genes
    
def find_chokepoint_reactions(map: str) -> list[str]:
    '''Find chokepoint reactions in the metabolic map

    Parameters
    ----------
    Metabolic map in SBML format (.xml)

    Returns
    -------
    Chokepoint reactions list
    '''
    
    model = CobraMetabolicModel(map)

    # Update flux bounds with FVA
    model.fva(update_flux = True)

    # Get chokepoints
    model.find_chokepoints(exclude_dead_reactions = True)
    chokepoints = model.chokepoints()

    chokepoint_list = []
    for i in range(len(chokepoints)):
        if chokepoints[i][0].id not in chokepoint_list:
            chokepoint_list.append(chokepoints[i][0].id)

    if initial_args.verbosity > 0:
        print(f"{clrs['g']}Finding chokepoint reactions...{clrs['n']}")
        model.print_chokepoints()

    # Remove duplicates
    chokepoint_list = list(set(chokepoint_list))

    # Get chokepoints
    return chokepoint_list

def find_essential_chokepoint_reactions(map: str) -> set[Gene]:
    '''Find essential genes reactions in the metabolic map

    Parameters
    ----------
    map (str)
        Metabolic map in SBML format (.xml)

    Returns
    -------
    set[cobra.Gene]
        Essencial genes reactions list
    '''
    
    model = CobraMetabolicModel(map)

    # Update flux bounds with FVA
    model.fva(update_flux = True)

    # Find chokepoints
    model.find_chokepoints(exclude_dead_reactions = True)
    chokepoints = model.chokepoints()

    # Find essential genes
    model.find_essential_genes_1()
    genes = list(model.essential_genes())
   
    # Get essential_chokepoint_reactions
    essential_CP = []
    for i in range(len(chokepoints)):
        for j in range(len(genes)):
            if chokepoints[i][0].gene_reaction_rule == genes[j].id:
                essential_CP.append(genes[j].id)

    # Remove duplicates
    essential_CP = list(set(essential_CP))

    if initial_args.verbosity > 0:
        print(f"{clrs['g']}Essential chokepoint genes...{clrs['n']}")
        print(essential_CP)
   
    return essential_CP
