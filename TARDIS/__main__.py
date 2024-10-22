#!/usr/bin/env python3

# Imports
###############################################################################
from TARDIS.Initialise import *

import TARDIS.FindTargets as FindTargets
#import TARDIS.HomologySearch as HomologySearch
import TARDIS.Toolbox as Toolbox

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

# Description
###############################################################################

# Functions
###############################################################################
   
# Main Function
###############################################################################
def main() -> None:
    '''Main function of the pipeline
    '''

    try:
        # Create the target directory    
        target_dir, target_subdir = Toolbox.create_dirs()

        # Get the metabolic map
        model_path = Toolbox.get_model(target_dir)

        # Create the metabolic model and update flux bounds with FVA (computing this only once makes the pipeline more efficient)
        model = FindTargets.initialize_model(model_path)

        # Find essential genes
        genes = FindTargets.find_essential_genes(model)

        # Save essential genes in a txt file with the same name as the input file
        _ = Toolbox.save_essential_genes(model_path, genes)
        
        # Find chokepoints
        chokepoints = FindTargets.find_chokepoint_reactions(model)

        # Find Essential chokepoint genes
        essential_CP = FindTargets.find_essential_chokepoint_reactions(model)

        # Create output dataframe
        _ = Toolbox.create_output_file(target_dir, genes, chokepoints, essential_CP)

        # If there are any essential chokepoint genes, retrieve their sequences
        if essential_CP:
            # Retrieve targets sequences
            _ = Toolbox.retrieve_targets(target_subdir, essential_CP)
        else:
            # Warn the user that no essential chokepoint genes were found
            print(f"{clrs['c']}INFO{clrs['n']}: No essential chokepoint genes found!")

        # Create database for specificity analysis
    except Exception as e:
        print(f"{clrs['r']}ERROR{clrs['n']}: {e}")
        exit(1)

# Execute
###############################################################################
if __name__ == "__main__":
    main()
