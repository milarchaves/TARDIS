#!/usr/bin/env python3

# Description
###############################################################################
'''
Module for performing the mapping analysis using CarveMe.

They are imported as:

import TARDIS.Mapping as Mapping
'''

# Imports
###############################################################################
from TARDIS.Initialise import *

import subprocess
import os

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

def create_map(input_path: str, output: str, template: str) -> str:
    '''Open input file, parse protein name and sequence to output dataframe and create the metabolic map using CarveMe (D. Machado et al, 2018. https://doi.org/10.1093/nar/gky537).

    Parameters
    ----------
    input_path (str)
        Input file in FASTA format (.faa)
    output (str)
        Output file in SBML format (.xml)
    template (str)
        Template to be used for the gram-positive bacteria in CarveMe reconstruction. Options: grampos or gramneg

    Returns
    -------
    str
        Path to metabolic map in SBML format (.xml)

    Raises
    ------
    ValueError
        No input file or input file in wrong format error
    subprocess.CalledProcessError
        Command failed with exit code error
    '''

    if template not in ["grampos", "gramneg"]:
        raise ValueError(f"Invalid template: {template}. Use 'grampos' or 'gramneg'.")

    # Get the filename from the input
    model = os.path.join(
        output,
        f"{os.path.splitext(os.path.basename(input_path))[0]}.xml"
    )

    # Create the command to run CarveMe
    command = [
        "carve", 
        input_path, 
        "--fbc2", 
        "--universe", 
        template, 
        "-o", 
        model,
        "--solver",
        initial_args.solver
    ]
        
    # Check if the input file is in FASTA format for DNA sequences
    if input_path.endswith(".fna"):
        # Add the DNA flag to the command
        command.append("--dna")

    # Check for verbosity
    if initial_args.verbosity > 0:
        print(f"{clrs['g']}Input file provided. Creating metabolic map...{clrs['n']}")
        print(f"{clrs['c']}Command{clrs['n']}: {' '.join(command)}")

    try:
        subprocess.run(
            command,
            check = True,
            text = True,
            capture_output = True
        )
    except subprocess.CalledProcessError as e:
        print(f"{clrs['r']}ERROR{clrs['n']}: Command failed with exit code {e.returncode}.")
        print(f"Error message: {e.stderr}")
        sys.exit(e.returncode)
    except Exception as e:
        print(f"{clrs['r']}ERROR{clrs['n']}: Unexpected error: {str(e)}")
        sys.exit(1)

    return model
