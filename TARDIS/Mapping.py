#!/usr/bin/env python3

# Description
###############################################################################
'''
Module for performing the mapping analysis using CarveMe.

They are imported as:

from TARDIS.Initialise import *
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
    '''Open input file, parse protein name and sequence to output dataframe and create the metabolic map using CarveMe (D. Machado et al, 2018. https://doi.org/10.1093/nar/gky537)

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
    No input file or input file in wrong format error
    '''

    # Get the filename from the input
    model = os.path.join(
        output,
        f"{os.path.splitext(os.path.basename(input_path))[0]}.xml"
    )

    if template in ["grampos", "gramneg"]:
        command = [
            "carve", 
            os.path.join(input_path), 
            "--fbc2", 
            "--universe", 
            template, 
            "-o", 
            model
        ]
    else:
        print(f"Unrecognized template. Expected either 'grampos' or 'gramneg' but got '{template}'.")
        exit(2)
    
    # Check if the input file is in FASTA format for DNA sequences
    if input_path.endswith(".fna"):
        # Add the DNA flag to the command
        command.append("--dna")

    # Check for verbosity
    if initial_args.verbosity > 0:
        print(f"{clrs['g']}Input file provided. Creating metabolic map...{clrs['n']}")
        print(f"Command: {' '.join(command)}")

    try:
        subprocess.run(
            command,
            check = True,
            text = True,
            capture_output = True
        )
    except subprocess.CalledProcessError as e:
        print(f"Command failed with exit code {e.returncode}.")
        print(f"Error message: {e.stderr}")
        print(f"Command: {' '.join(command)}")
        exit(e.returncode)
    except Exception as e:
        print(f"Unexpected error: {str(e)}")
        print(f"Command: {' '.join(command)}")
        exit(1)

    return model