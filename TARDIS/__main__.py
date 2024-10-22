#!/usr/bin/env python3

# Imports
###############################################################################
from TARDIS.Initialise import *

import TARDIS.FindTargets as FindTargets
#import TARDIS.HomologySearch as HomologySearch
import TARDIS.Mapping as Mapping

import os
import pandas as pd
from Bio import SeqIO

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

# Global variables
###############################################################################
global output

# Functions
###############################################################################

# Main Function
###############################################################################
def main() -> None:
    '''Main function of the pipeline
    '''

    # Check if the output directory exists
    if not os.path.isdir(initial_args.output):
        try:
            os.mkdir(initial_args.output)
        except OSError as e:
            print(f"Creation of the directory {initial_args.output} failed. Error: {e}")
            exit(3)
        except Exception as e:
            print(f"Unknown error occurred: {e}")
            exit(4)

    # Get the target directory name
    target_dir = os.path.join(
        initial_args.output,
        os.path.splitext(os.path.basename(initial_args.input_file))[0]
    )
    
    # Check if the target directory exists
    if not os.path.isdir(target_dir):
        # Create a directory to save the targets sequences
        os.mkdir(target_dir)

    # Check if the input file and the template are provided
    if initial_args.input_file and initial_args.template:
        # Create metabolic map
        model = Mapping.create_map(initial_args.input_file, initial_args.output, initial_args.template)
    elif initial_args.metabolic_map:
        # Use the provided metabolic map
        model = initial_args.metabolic_map
    else:
        print("Please, provide a metabolic map in SBML format or a proteome file indicating the template (grampos or gramneg) to be used to build the pyth. For more information, use -h or --help")
        exit(5)

    # Find essential genes
    genes = FindTargets.find_essential_genes(model)

    # Save essential genes in a txt file with the same name as the input file
    with open(model.replace(".xml", ".txt"), 'w') as f:
        for item in genes:
            f.write(f"{item}\n")
    
    # Find chokepoints
    chokepoints = FindTargets.find_chokepoint_reactions(model)

    # Find Essential chokepoint genes
    essential_CP = FindTargets.find_essential_chokepoint_reactions(model) 
    
    genes_df = pd.Series(list(genes))
    chokepoints_df = pd.Series(list(chokepoints))
    essential_CP_df = pd.Series(list(essential_CP))

    frame = {
        "Essential genes": genes_df, 
        "Chokepoint Reactions": chokepoints_df, 
        "Essential chokepoint genes": essential_CP_df
    }

    output = pd.DataFrame(frame)

    # Save output dataframe
    output.to_csv(
        os.path.join([initial_args.output,  "output.csv"]),
        sep=',',
        index = False
    )
    
    # Retrieve targets sequences
    if initial_args.input_file:
        with open(initial_args.input_file, 'r'):
            seqs = list(SeqIO.parse(initial_args.input_file, "fasta"))

        for seq in seqs:
            if seq.id.replace('|', '_') in essential_CP:
                try:
                    name = os.path.join("targets", seq.id.split('|')[1] + ".fasta")
                    
                    with open(name, 'a') as f:
                        f.write(f">{seq.id}\n{str(seq.seq)}")

                    if initial_args.verbosity > 0:
                        print(f"{clrs['g']}Targets sequences...{clrs['n']}")
                        print(seq.id)
                        print(seq.seq)
                        print('\n')            
                except:
                    pass
    
    # Create database for specificity analysis

# Execute
###############################################################################
if __name__ == "__main__":
    main()
