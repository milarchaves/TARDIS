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


def create_tardis_conf() -> None:
    '''Creates the 'TARDIS.conf' file.

    Parameters
    ----------
    None

    Returns
    -------
    None

    Raises
    ------
    None
    '''
    conf_file = "TARDIS.conf"
    with open(conf_file, 'w') as cf:
        cf.write(tw.dedent("""
        # CarveMe Executable
        carveme_exe = "/home/camila/anaconda3/envs/tardis/bin/carve"

        # FindCP Executable
        findcp_exe = """))

    print(clrs['g']+'Configuration file created!'+clrs['n']+' Please'+clrs['y']+' EDIT ITS CONTENTS '+clrs['n']+'to match your environment and run TARDIS again.')

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


    parser.add_argument('-f', '--file',
                        dest='input_file',
                        type=str,
                        metavar='',
                        help='File containing the paths of bacterial genome .fna or proteome .faa')

    parser.add_argument('-o', '--output',dest='output',
                        type=str,
                        metavar='',
                        help='Defines the output path')

    parser.add_argument('-v', '--verbose',
                        dest='verbosity',
                        action='count',
                        default=0,
                        help='Controls verbosity')

    parser.add_argument('--conf',
                        dest='config_file',
                        type=str,
                        metavar='',
                        help='Configuration file containing external executable paths')

    parser.add_argument('-vn', '--version', action='version',
                    version=f'%(prog)s {TARDISVersion}')

    initial_args = parser.parse_args()

    return initial_args

initial_args = argument_parsing()

print(description)





# Criar classe proteína com nome, gene e sequência 