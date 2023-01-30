#!/usr/lib/python3

# Imports
###############################################################################
from TARDIS.Initialise import *
import subprocess
import sys
import os
import pandas as pd

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

# Classes
###############################################################################

# Functions
###############################################################################

# Main Function
###############################################################################
def main():
    output = pd.DataFrame(columns=['Protein', 'Sequence', 'Essencial gene', 'Chokepoint', 'Homologue'])
    # Open input file and create the metabolic map
    if initial_args.input_file is None:
        print(clrs['r']+'ERROR: '+clrs['n']+'No input file was provided. Please use the'+clrs['y']+' -f '+clrs['n']+'flag to provide a file with the paths of the genomes/proteomes.')
    print(initial_args.input_file)
    if initial_args.input_file.endswith('.faa'):
        if initial_args.verbosity > 0:
            print(clrs['g']+'Input file provided. Creating metabolic map...'+clrs['n'])
        print(" ".join(["/home/camila/anaconda3/envs/tardis/bin/carve", initial_args.input_file]))
        subprocess.run(["carve", initial_args.input_file], capture_output=True)
        with open(initial_args.input_file, 'r') as f:
            for line in f:
                sequence = ''
                if line[0]=='>':
                    name = line.split()[0].replace('>', '').replace(' ','_')
                else:
                    while line[0]!='>':
                        sequence += line
                        line = f.readline()
                        if not line:
                            break
                output = pd.concat([output, pd.DataFrame([[name, sequence, '', '', '']], columns=['Protein', 'Sequence', 'Essencial gene', 'Chokepoint', 'Homologue'])], ignore_index=True)
                
    output.to_csv('output.csv', index=False)
# Execute
###############################################################################
if __name__ == "__main__":
    main()
