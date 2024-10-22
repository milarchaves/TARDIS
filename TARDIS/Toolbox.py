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
        os.makedirs(dirname, exist_ok = True)
        
        if initial_args.verbosity > 0:
            print(f"{clrs['g']}Directory created{clrs['n']}: {dirname}")
    # Catch the exception if the directory cannot be created
    except OSError as e:
        print(f"Creation of the directory {dirname} failed. Error: {e}")
        exit(3)

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

    # Get the target directory name
    target_dir = os.path.join(
        initial_args.output,
        os.path.splitext(os.path.basename(initial_args.input_file))[0]
    )

    # Parameterize the target directory
    target_subdir = os.path.join(target_dir, "targets")

    # Create the output directory
    create_dir(initial_args.output)
    # Create the target directory
    create_dir(target_dir)
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

    frame = {
        "Essential genes": pd.Series(list(genes)), 
        "Chokepoint Reactions": pd.Series(list(chokepoints)), 
        "Essential chokepoint genes": pd.Series(list(essential_CP))
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

    # Check if the input file is provided
    if not initial_args.input_file:
        return
    
    try:
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
        # Parameterize the gene id
        gene_id = seq.id.replace('|', '_')

        # If the sequence is in the essential chokepoint genes list
        if gene_id in essential_CP: # TODO: make this more robust
            write_sequence_to_file(target_subdir, seq)
    
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

def write_sequence_to_file(target_subdir: str, seq) -> None:
    '''Write a sequence to a FASTA file.

    Parameters
    ----------
    target_subdir : str
        Path to the target subdirectory

    seq : SeqRecord
        Sequence record
    '''

    try:
        # Get the file name
        file_name = f"{seq.id.split('|')[1]}.fasta" if '|' in seq.id else f"{seq.id}.fasta"

        # Get the file path
        file_path = os.path.join(target_subdir, file_name)

        with open(file_path, 'a') as f:
            f.write(f">{seq.id}\n{str(seq.seq)}")

        if initial_args.verbosity > 0:
            print(f"{clrs['g']}Targets sequences...{clrs['n']}")
            print(seq.id)
            print(seq.seq)
            print('\n')            
    except:
        print(f"{clrs['r']}ERROR{clrs['n']}: An error occurred while writing the sequence '{seq.id}' to the file.")
