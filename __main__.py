#!/usr/lib/python3

# Imports
###############################################################################
from TARDIS.Initialise import *
from FindTargets import *
import subprocess
import sys
import os
import pandas as pd

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

# Description
###############################################################################

# Global variables
###############################################################################
global output
# Functions
###############################################################################

# Main Function
###############################################################################
def main():
    # Create dataframe output
    output =pd.DataFrame(columns=['Essential genes', 'Chokepoint Reactions', 'Essential chokepoint genes'])
    #create metabolic map
    input = initial_args.input_file
    
    map = create_map(input)
    #find targets
    genes = find_essential_genes(map)
    
    #find chokepoints
    chokepoints = find_chokepoint_reactions(map)

    #find Essential chokepoint genes
    essential_CP =find_essential_chokepoint_reactions(map) 
    
    genes = list(genes)
    chokepoints = list(chokepoints)
    essential_CP = list(essential_CP)

    frame = {'Essential genes': genes, 'Chokepoint Reactions': chokepoints, 'Essential chokepoint genes': essential_CP}

    output =pd.DataFrame(frame)
    
    #save output dataframe
    output.to_csv('output.csv', sep='\t', index=False)

# Execute
###############################################################################
if __name__ == "__main__":
    main()
