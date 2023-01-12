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

# Description
###############################################################################
description = tw.dedent("""
    \033[1;93m+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
    +-+-+-+-+-+-+-+-+-+-+-+-+-+  \033[1;95m   _____  _    ____  ____ ___ ____  
|_   _|/ \  |  _ \|  _ \_ _/ ___| 
  | | / _ \ | |_) | | | | |\___ \ 
  | |/ ___ \|  _ <| |_| | | ___) |
  |_/_/   \_\_| \_\____/___|____/ 
                                  
    \033[1;93m+-+-+-+-+-+-+-+-+-+-+-+-+-+
    +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
    \033[1;95m
                  “There’s always something to look at if you open your eyes!”
    \033[1;93m
                                 - The Fifth Doctor
    \033[1;0m
      Copyright (C) 2022  Chaves, C; Torres, P.H.M.

    \033[1;93m
                      [Federal University of Rio de Janeiro]

    \033[1;0m
          This program comes with ABSOLUTELY NO WARRANTY

      TARDIS uses metabolic network moddeling to identify potencial targets for drug design against infectious diseases. 
     \033[1;0m+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
          """)

epilogue = tw.dedent("""
    TARDIS generates a summary file, whose columns are ordered as follows:

    1 Protein Name  
    2 Protein Sequence
    3 Human Homologue
    4 Identity
    """)

# functions
###############################################################################


def create_tardis_conf():
    conf_file = "tardis.cfg"
    with open(conf_file, 'w') as cf:
        cf.write(tw.dedent("""
        # CarveMe Executable
        carveme_exe = 

        # FindCP Executable
        findcp_exe = """))

    print(clrs['g']+'Configuration file created!'+clrs['n']+' Please'+clrs['y']+' EDIT ITS CONTENTS '+clrs['n']+'to match your environment and run TARDIS again.')

# Define Global Variables
###############################################################################

# Parse command line arguments
###############################################################################
def argument_parsing():
    parser = argparse.ArgumentParser(prog='TARDIS',
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=description,
                                     epilog=epilogue)

    parser.add_argument('--version', action='version',
                    version='%(prog)s 1.0')

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

    initial_args = parser.parse_args()

    return initial_args

initial_args = argument_parsing()




# Criar classe proteína com nome, gene e sequência 