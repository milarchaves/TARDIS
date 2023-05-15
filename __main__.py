#!/usr/lib/python3

# Imports
###############################################################################
from TARDIS.Initialise import *
from FindTargets import *
import os
import pandas as pd
from Bio import SeqIO

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
    #create metabolic map
    input = initial_args.input_file
    
    map = create_map(input)

    #find essential genes
    genes = find_essential_genes(map)
    #save essential genes in a txt file with the same name as the input file
    with open(input.replace(".fasta", ".txt"), 'w') as f:
        for item in genes:
            f.write("%s\n" % item)
    
    #find chokepoints
    chokepoints = find_chokepoint_reactions(map)

    #find Essential chokepoint genes
    essential_CP =find_essential_chokepoint_reactions(map) 
    
    genes_df = pd.Series(list(genes))
    chokepoints_df = pd.Series(list(chokepoints))
    essential_CP_df = pd.Series(list(essential_CP))

    frame = {'Essential genes': genes_df, 'Chokepoint Reactions': chokepoints_df, 'Essential chokepoint genes': essential_CP_df}

    output = pd.DataFrame(frame)
    #save output dataframe
    output.to_csv('output.csv', sep='\t', index=False)
    
    #Retrieve targets sequences
    with open(input, 'r'):
        seqs = list(SeqIO.parse(input, "fasta"))
    
    #create a directory to save the targets sequences
    namedir = input.replace(".fasta", "")
    if not os.path.isdir(namedir):
        os.mkdir(namedir)
    for seq in seqs:
        if seq.id.replace('|', '_') in essential_CP:
            try:
                name = os.path.join('targets', seq.id.split('|')[1] + '.fasta')
                with open(name, 'a') as f:
                    f.write('>' + seq.id + '\n' + str(seq.seq))
                if initial_args.verbosity > 0:
                    print(clrs['g']+'Targets sequences...'+clrs['n'])
                    print(seq.id)
                    print(seq.seq)
                    print('\n')            
            except:
                pass
  
# Execute
###############################################################################
if __name__ == "__main__":
    main()
