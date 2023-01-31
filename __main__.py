#!/usr/lib/python3

# Imports
###############################################################################
from TARDIS.Initialise import *
from FindTargets import *
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

# Global variables
###############################################################################
global output
# Functions
###############################################################################

# Main Function
###############################################################################
def main():
    # Create dataframe output
    output =pd.DataFrame(columns=['Protein', 'Sequence', 'Essential gene', 'Chokepoint', 'Homologue'])
    #create metabolic map
    input = initial_args.input_file
    
    map = create_map(input)
    #find targets
    genes = find_essential_genes(map)
    #find chokepoints
    chokepoints = find_chokepoint_reactions(map)

    
    
    genes = list(genes)
    for i in range(len(genes)):
        output = pd.concat([output, pd.DataFrame([['', '', genes[i], '', '']], columns=['Protein', 'Sequence', 'Essencial gene', 'Chokepoint', 'Homologue'])], ignore_index=True)

    #update the output dataframe
    with open(input, 'r') as f:
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
                if name in output['Essencial gene']:
                    output = pd.concat([output, pd.DataFrame([[name, sequence, '', '', '']], columns=['Protein', 'Sequence', 'Essencial gene', 'Chokepoint', 'Homologue'])], ignore_index=True)

    print(output)            
    '''
    for gene in model.genes:
        print (gene)
        if gene.id in essential_genes:
            output.loc[output['Protein'] == gene.id, 'Essential gene'] = 'Yes'
        if gene.id in chokepoints:
            output.loc[output['Protein'] == gene.id, 'Chokepoint'] = 'Yes'
    '''
    #save output dataframe
    #output.to_csv(output, sep='\t', index=False)

# Execute
###############################################################################
if __name__ == "__main__":
    main()
