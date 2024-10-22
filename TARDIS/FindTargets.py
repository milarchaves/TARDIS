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
from functools import lru_cache

from typing import Union

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

def compute_if_needed(model: CobraMetabolicModel, methods: list[str]) -> None:
    '''Ensure essential genes or chokepoints are computed.

    Parameters
    ----------
    model : CobraMetabolicModel
        CobraMetabolicModel object.
    methods : list[str]
        List of methods to compute.
    
    Raises
    ------
    ValueError
        If methods is not a list.
    '''

    # Check if methods is a list
    if not isinstance(methods, list):
        raise ValueError("Provide a list of methods to compute.")

    # For each method, check if it is already computed
    for method in methods:
        if method == 'essential_genes':
            # Try to retrieve essential genes, if it fails, compute them
            try:
                _ = model.essential_genes()
            except:
                model.find_essential_genes_1()
        elif method == 'chokepoints':
            # Try to retrieve chokepoints, if it fails, compute them
            try:
                _ = model.chokepoints()
            except:
                model.find_chokepoints(exclude_dead_reactions=True)

def find_essential_genes(model: Union[str, CobraMetabolicModel]) -> set[Gene]:
    '''Find essential genes in the metabolic map

    Parameters
    ----------
    model (CobraMetabolicModel | str)
        Path to the metabolic map in SBML format (.xml) or CobraMetabolicModel object.

    Returns
    -------
    set[cobra.Gene]
        Essencial genes list
    
    Raises
    ------
    ValueError
        If model is not a path to a metabolic map in SBML format (.xml) or a CobraMetabolicModel object.
    '''
    
    # Load or validate the metabolic model
    model = get_model(model)

    # Check if the essential genes are already computed
    compute_if_needed(model, ["essential_genes"])
    
    if initial_args.verbosity > 0:
        print(f"{clrs['g']}Finding essential genes...{clrs['n']}")
        model.print_essential_genes()

    return model.essential_genes()
    
def find_chokepoint_reactions(model: Union[str, CobraMetabolicModel]) -> list[str]:
    '''Find chokepoint reactions in the metabolic map

    Parameters
    ----------
    model (CobraMetabolicModel | str)
        Path to the metabolic map in SBML format (.xml) or CobraMetabolicModel object.

    Returns
    -------
    list[str]
        Chokepoint reactions list.
    '''
    
    # Load or validate the metabolic model
    model = get_model(model)

    # Check if the chokepoints are already computed
    compute_if_needed(model, ["chokepoints"])

    if initial_args.verbosity > 0:
        print(f"{clrs['g']}Finding chokepoint reactions...{clrs['n']}")
        model.print_chokepoints()
        
    # Get chokepoints using set comprehension to remove duplicates
    return list({cp[0].id for cp in model.chokepoints()})

def find_essential_chokepoint_reactions(model: Union[str, CobraMetabolicModel]) -> list[str]:
    '''Find essential genes reactions in the metabolic map

    Parameters
    ----------
    model (CobraMetabolicModel | str)
        Path to the metabolic map in SBML format (.xml) or CobraMetabolicModel object.

    Returns
    -------
    list[str]
        Essencial genes reactions ids list
    '''
    
    # Load or validate the metabolic model
    model = get_model(model)

    # Check if the chokepoints are already computed
    compute_if_needed(model, ["essential_genes", "chokepoints"])

    # Get chokepoints using set comprehension to remove duplicates
    chokepoints = {cp[0].gene_reaction_rule for cp in model.chokepoints()}

    # Get essential genes using set comprehension to remove duplicates
    essential_genes = {gene.id for gene in model.essential_genes()}
   
    # Get essential_chokepoint_reactions
    essential_CP = list(chokepoints.intersection(essential_genes))

    if initial_args.verbosity > 0:
        print(f"{clrs['g']}Essential chokepoint genes...{clrs['n']}")
        print(essential_CP)
   
    return essential_CP

def get_model(model: Union[str, CobraMetabolicModel]) -> CobraMetabolicModel:
    '''Load or validate the metabolic model.

    Parameters
    ----------
    model : str | CobraMetabolicModel
        Path to the metabolic map in SBML format (.xml) or CobraMetabolicModel object.

    Returns
    -------
    CobraMetabolicModel
        Initialized metabolic model with FVA applied.
    '''

    # Check if the model is a path to a metabolic map in SBML format (.xml)
    if isinstance(model, str):
        print(f"{clrs['y']}WARNING{clrs['n']}: Use a CobraMetabolicModel object for efficiency.")
        return initialize_model(model)
    
    # Check if the model is a CobraMetabolicModel object
    if not isinstance(model, CobraMetabolicModel):
        raise ValueError("Provide a valid SBML path or a CobraMetabolicModel object.")
    
    return model

# Cache the model initialization to avoid recomputing it
@lru_cache
def initialize_model(map: str) -> CobraMetabolicModel:
    '''Initialize the metabolic model with FVA applied.

    Parameters
    ----------
    map : str
        Path to the metabolic map in SBML format (.xml)

    Returns
    -------
    CobraMetabolicModel
        Initialized metabolic model with FVA applied.
    '''

    # Create the metabolic model
    model = CobraMetabolicModel(map)

    # Update flux bounds with FVA
    model.fva(update_flux = True)

    return model
