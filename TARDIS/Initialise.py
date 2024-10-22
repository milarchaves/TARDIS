#!/usr/bin/env python3

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
import os

from reframed import set_default_solver

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
def argument_parsing() -> argparse.Namespace:
    '''Parse command line arguments

    Parameters
    ----------
    None

    Returns
    -------
    argparse.Namespace
        Parsed arguments
    '''

    parser = argparse.ArgumentParser(prog = "TARDIS",
                                     formatter_class = argparse.RawDescriptionHelpFormatter,
                                     description = description,
                                     epilog = epilogue)

    parser.add_argument("-i", "--input",
                        dest = "input_file",
                        type = str,
                        help = "Path of bacterial genome (.fna) or proteome (.faa) file.")

    parser.add_argument("-m", "--metabolic-map",
                        dest = "metabolic_map",
                        type = str,
                        help = "Path of metabolic map in SBML format (.xml).")

    parser.add_argument("-t", "--template",
                        dest = "template",
                        type = str,
                        help = "Indicate the template to be used for the gram-positive bacteria in CarveMe reconstruction. grampos for gram-positive bacteria and gramneg for gram-negative bacteria.")

    parser.add_argument("-s", "--solver",
                        dest = "solver",
                        type = str,
                        default = "cplex",
                        help = "Select MILP solver. Available options: cplex [default], gurobi, scip.")

    parser.add_argument("-o", "--output",
                        dest = "output",
                        type = str,
                        help = "Defines the output path.")

    parser.add_argument("-v", "--verbose",
                        dest = "verbosity",
                        action = "count",
                        default = 0,
                        help = "Controls verbosity.")

    parser.add_argument("-vn", "--version", 
                        action = "version",
                        version = f"%(prog)s {TARDISVersion}")

    initial_args = parser.parse_args()

    return initial_args

# Splash, version & clear tmp
###############################################################################
TARDISVersion = "1.0"

description = tw.dedent("""\033[1;96m
    +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-++-+-+-+-+
    +-+-+-+-+-+-+-+-+-+-\033[1;95m  _____  _    ____  ____ ___ ____   \033[1;96m-+-+-+-+-+-+-+-+-+-+  
    +-+-+-+-+-+-+-+-+-+-\033[1;95m |_   _|/ \  |  _ \|  _ \_ _/ ___|  \033[1;96m-+-+-+-+-+-+-+-+-+-+ 
    +-+-+-+-+-+-+-+-+-+-\033[1;95m   | | / _ \ | |_) | | | | |\___    \033[1;96m-+-+-+-+-+-+-+-+-+-+ 
    +-+-+-+-+-+-+-+-+-+-\033[1;95m   | |/ ___ \|  _ <| |_| | | ___) | \033[1;96m-+-+-+-+-+-+-+-+-+-+
    +-+-+-+-+-+-+-+-+-+-\033[1;95m   |_/_/   \_\_| \_\____/___|____/  \033[1;96m-+-+-+-+-+-+-+-+-+-+
    +-+-+-+-+-+-+-+-+-+-\033[1;95m                                    \033[1;96m-+-+-+-+-+-+-+-+-+-+
    +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-++-+-+-+-+
\033[1;0m
      “There’s always something to look at if you open your eyes!”
\033[1;0m
                                                - The Fifth Doctor
\033[1;0m
      Copyright (C) 2022  Chaves, C; Rossi, A.D; Torres, P.H.M.
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
    "n": "\033[1;0m"    # default
}

# Parse the arguments
initial_args = argument_parsing()

# Print splash
print(description)

# Solver checking
try:
    # Try to set the default solver
    set_default_solver(initial_args.solver)
except Exception as e1:
    # Check if the chosen solver is not scip
    if initial_args.solver != "scip":
        # Warn the user that the solver will be changed to SCIP (and it will take longer to run)
        print(f"{clrs['y']}WARNING{clrs['n']}: The solver {initial_args.solver} is not available. Using SCIP instead. The analysis may take longer to complete.")
        
        try:
            # Set the default solver to SCIP
            initial_args.solver = "scip"

            # Try to set the default solver to SCIP
            set_default_solver(initial_args.solver)
        except Exception as e2:
            # If it fails, print an error message and exit
            print(f"{clrs['r']}ERROR{clrs['n']}: {e2}")
            exit(1)
    else:
        print(f"{clrs['r']}ERROR{clrs['n']}: {e1}")
        exit(1)
