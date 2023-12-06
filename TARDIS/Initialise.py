#!/usr/lib/python3

# Description
###############################################################################
'''
This file contains variables and functions that are used to initialise TARDIS.\n
All scripts that use TARDIS must import this file.

They are imported as:

from TARDIS.Initialise import *
'''

# Imports
###############################################################################
import argparse
import textwrap as tw
import subprocess
import os

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
# Splash, version & clear tmp
###############################################################################
TARDISVersion = "1.0"

description = tw.dedent("""\033[1;96m
    +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-++-+-+-+-+
    +-+-+-+-+-+-+-+-+-+-\033[1;95m  _____  _    ____  ____ ___ ____   \033[1;96m-+-+-+-+-+-+-+-+-+-+  
    +-+-+-+-+-+-+-+-+-+-\033[1;95m |_   _|/ \  |  _ \|  _ \_ _/ ___|  \033[1;96m-+-+-+-+-+-+-+-+-+-+ 
    +-+-+-+-+-+-+-+-+-+-\033[1;95m   | | / _ \ | |_) | | | | |\___    \033[1;96m-+-+-+- +-+-+-+-+-+-+ 
    +-+-+-+-+-+-+-+-+-+-\033[1;95m   | |/ ___ \|  _ <| |_| | | ___) | \033[1;96m-+-+-+-+-+-+-+-+-+-+
    +-+-+-+-+-+-+-+-+-+-\033[1;95m   |_/_/   \_\_| \_\____/___|____/  \033[1;96m-+-+-+-+-+-+-+-+-+-+
    +-+-+-+-+-+-+-+-+-+-\033[1;95m                                    \033[1;96m-+-+-+-+-+-+-+-+-+-+
    +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-++-+-+-+-+
\033[1;0m
      “There’s always something to look at if you open your eyes!”
\033[1;0m
                                                - The Fifth Doctor
\033[1;0m
      Copyright (C) 2022  Chaves, C; Torres, P.H.M.
\033[1;91m
                  [The Federal University of Rio de Janeiro]
\033[1;0m
                This program comes with ABSOLUTELY NO WARRANTY
\033[1;0m
      TARDIS uses genome-scale metabolic network moddeling to find new targets for antimicrobial drug-design
\033[1;0m""")

epilogue = tw.dedent("""
    TARDIS generates a summary file, whose columns are ordered as follows:

    1 Protein Name  
    2 Protein Sequence
    3 Human Homologue
    4 Identity
    """)

# functions
###############################################################################

def create_map (input, template) -> str:
    '''Open input file, parse protein name and sequence to output dataframe and create the metabolic map

    Parameters
    ----------
    Input file in FASTA format (.faa)

    Returns
    -------
    Metabolic map in SBML format (.xml)

    Raises
    ------
    No input file or input file in wrong format error
    '''
    if input.endswith('.faa') or input.endswith('.fasta'):
        if initial_args.verbosity > 0:
            print(clrs['g']+'Input file provided. Creating metabolic map...'+clrs['n'])
        if template == 'grampos':
            subprocess.run(["carve", os.path.join(input), "--fbc2", "--universe", "grampos"])
        elif template == 'gramneg':
            subprocess.run(["carve", os.path.join(input), "--fbc2", "--universe", "gramneg"])
        if input.endswith('.faa'):
            model = str(input).replace('.faa', '.xml')
        else:
            model = str(input).replace('.fasta', '.xml')
    elif input.endswith('.fna'):
        if initial_args.verbosity > 0:
            print(clrs['g']+'Input file provided. Creating metabolic map...'+clrs['n'])
        if template == 'grampos':
            subprocess.run(["carve", "--dna", os.path.join(input), "--fbc2", "--universe", "grampos"])
        elif template == 'gramneg':
            subprocess.run(["carve", "--dna", os.path.join(input), "--fbc2", "--universe", "gramneg"])
        model = str(input).replace('.fna', '.xml')

    else:
        print(clrs['r']+'ERROR: '+clrs['n']+'The input file must be in FASTA format')
    return model


# Define Global Variables
###############################################################################

# Dictionary for the output colors
clrs = {
    "r": "\033[1;91m",  # red
    "g": "\033[1;92m",  # green
    "y": "\033[1;93m",  # yellow
    "b": "\033[1;94m",  # blue
    "p": "\033[1;95m",  # purple
    "c": "\033[1;96m",  # cyan
    "n": "\033[1;0m"   # default
    }

# Parse command line arguments
###############################################################################
def argument_parsing():
    parser = argparse.ArgumentParser(prog='TARDIS',
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=description,
                                     epilog=epilogue)


    parser.add_argument('-i', '--input',
                        dest='input_file',
                        type=str,
                        help='Path of bacterial genome (.fna) or proteome (.faa) file')
    
    parser.add_argument('-m', '--metabolic-map',
                        dest='metabolic_map',
                        type=str,
                        help='Path of metabolic map in SBML format (.xml)')
    
    parser.add_argument('-t', '--template',
                        dest='template',
                        type=str,
                        help='Indicate the template to be used for the gram-positive bacteria in CarveMe reconstruction. grampos for gram-positive bacteria and gramneg for gram-negative bacteria')

    parser.add_argument('-o', '--output',
                        dest='output',
                        type=str,
                        help='Defines the output path')
    
    parser.add_argument('-d', '--database',
                        dest='database',
                        type=str,
                        help='Creates database for specificity analysis from a proteome file in fasta format')
    
    parser.add_argument('-dn', '--database-name',
                        dest='database_name',
                        type=str,
                        help='Names the database for specificity analysis')
    
    parser.add_argument('-e', '--evalue',
                    dest='evalue',
                    type=float,
                    default=0.04,
                    help='Defines the e-value threshold for psi-blast (default: 0.04)')

    parser.add_argument('-v', '--verbose',
                        dest='verbosity',
                        action='count',
                        default=0,
                        help='Controls verbosity')

    parser.add_argument('-vn', '--version', 
                        action='version',
                        version=f'%(prog)s {TARDISVersion}')

    initial_args = parser.parse_args()

    return initial_args

initial_args = argument_parsing()

print(description)





# Criar classe proteína com nome, gene e sequência 