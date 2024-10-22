#!/usr/bin/env python3

# Description
###############################################################################
'''
Various functions that are used in the TARDIS scripts.

They are imported as:

import TARDIS.Toolbox as Toolbox
'''

# Imports
###############################################################################
from TARDIS.Initialise import *

import TARDIS.Mapping as Mapping

from Bio import SeqIO
from cobra import Gene

import os
import pandas as pd

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

# Functions
###############################################################################
def create_dir(dirname: str) -> None:
    ''' Create a directory

    Parameters
    ----------
    dirname : str
        The name of the directory to be created
    '''

    try:
        # If the directory does not exist, create it
        if not os.path.exists(dirname):
            os.mkdir(dirname)
        
        if initial_args.verbosity > 0:
            print(f"{clrs['g']}Directory created{clrs['n']}: {dirname}")
    # Catch the exception if the directory cannot be created
    except OSError as e:
        print(f"Creation of the directory {dirname} failed. Error: {e}")
        exit(3)
    # Catch any other exception
    except Exception as e:
        print(f"Unknown error occurred: {e}")
        exit(4)
    
    return None

def create_dirs() -> tuple[str, str]:
    ''' Create all the necessary directories.

    Parameters
    ----------
    None

    Returns
    -------
    tuple(str, str)
        The target directory and its subdirectory
    '''

    # Create the output directory
    create_dir(initial_args.output)

    # Get the target directory name
    target_dir = os.path.join(
        initial_args.output,
        os.path.splitext(os.path.basename(initial_args.input_file))[0]
    )

    # Create the target directory
    create_dir(target_dir)

    # Parameterize the target directory
    target_subdir = os.path.join(target_dir, "targets")

    # Create the target subdirectory
    create_dir(target_subdir)
    
    return target_dir, target_subdir

def create_output_file(target_dir: str, genes: set[Gene], chokepoints: list[str], essential_CP: set[Gene]) -> None:
    ''' Create the output file

    Parameters
    ----------
    target_dir : str
        Path to the target directory
    genes : set[Gene]
        Essential genes list
    chokepoints : list[str]
        Chokepoint reactions list
    essential_CP : set[Gene]
        Essential chokepoint genes list
    '''

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
        os.path.join(target_dir, "output.csv"),
        sep=',',
        index = False
    )

    return None

def get_model(target_dir: str) -> str:
    ''' Get the metabolic map

    Parameters
    ----------
    target_dir : str
        Path to the target directory

    Returns
    -------
    str
        Path to the metabolic map in SBML format (.xml)
    '''

    # Check if the input file and the template are provided
    if initial_args.input_file and initial_args.template:
        # Create metabolic map
        return Mapping.create_map(initial_args.input_file, target_dir, initial_args.template)
    elif initial_args.metabolic_map:
        # Use the provided metabolic map
        return initial_args.metabolic_map
    else:
        print("Please, provide a metabolic map in SBML format or a proteome file indicating the template (grampos or gramneg) to be used to build the pyth. For more information, use -h or --help")
        exit(5)

def retrieve_targets(target_subdir: str, essential_CP: set[Gene]) -> None:
    ''' Retrieve targets sequences

    Parameters
    ----------
    target_subdir : str
        Path to the target subdirectory
    essential_CP : set[Gene]
        Essential chokepoint genes list
    '''

    # If the input file has been provided
    if initial_args.input_file:
        try:
            # Open the input file
            with open(initial_args.input_file, 'r'):
                # Parse the sequences
                seqs = list(SeqIO.parse(initial_args.input_file, "fasta"))
        except FileNotFoundError as e:
            print(f"{clrs['r']}ERROR{clrs['n']}: The input file was not found.")
            print(f"Error message: {e}")
            exit(6)
        except Exception as e:
            print(f"{clrs['r']}ERROR{clrs['n']}: An unknown error occurred.")
            print(f"Error message: {e}")
            exit(7)

        # For each sequence in the input file
        for seq in seqs:
            # If the sequence is in the essential chokepoint genes list
            if seq.id.replace('|', '_') in essential_CP: # TODO: make this more robust
                try:
                    if '|' in seq.id:
                        name = os.path.join(target_subdir, seq.id.split('|')[1] + ".fasta")
                    else:
                        name = os.path.join(target_subdir, seq.id + ".fasta")
                    
                    with open(name, 'a') as f:
                        f.write(f">{seq.id}\n{str(seq.seq)}")

                    if initial_args.verbosity > 0:
                        print(f"{clrs['g']}Targets sequences...{clrs['n']}")
                        print(seq.id)
                        print(seq.seq)
                        print('\n')            
                except:
                    pass

    return None

def save_essential_genes(model: str, genes: set[Gene]) -> None:
    ''' Save essential genes in a txt file with the same name as the input file

    Parameters
    ----------
    model : str
        Path to the metabolic map in SBML format (.xml)
    genes : set[Gene]
        Essential genes list
    '''

    # Save essential genes in a txt file with the same name as the input file
    with open(model.replace(".xml", ".txt"), 'w') as f:
        for item in genes:
            f.write(f"{item}\n")
    
    return None
