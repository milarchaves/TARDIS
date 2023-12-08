#!/usr/lib/python3

# Description
###############################################################################
'''
This file implements the protocol to look for specificity in the targets found. 

1. Given a FASTA file with protein sequences, creates a local database
2. Performs a Psi-BLAST search against the database
3. Returns targets with similarity below a given e-value threshold (default: 0.04)
  
'''

# Imports
###############################################################################
from Initialise import *
import os
import subprocess
from Bio.Blast import NCBIXML

# Functions
###############################################################################
def make_database(fasta_file, database_name):
    # make a local database from the fasta file
    make_db_cmd = ['makeblastdb',
                   '-in', fasta_file,
                   '-dbtype', 'prot',
                   '-out', database_name,
                   '-title', database_name,
                   '-parse_seqids']
    subprocess.run(make_db_cmd, shell=True)

def blast_check(database, targets, evalue):

    # set the output file name
    output_file = "specificity_results.xml"

    # perform psi-blast on all fasta files in the directory
    for out_dir in os.listdir(targets):
        if out_dir.endswith('.fasta'):
            print("Blasting " + out_dir)
            query_file = os.path.join(targets, out_dir)
            blast_cmd = ['psiblast',
                 '-query', query_file,
                 '-db', database,
                 '-num_iterations', "3",
                 '-out', output_file,
                 '-outfmt', '5']
            subprocess.run(blast_cmd)
            #print(" ".join(blast_cmd))

    # Iterate through the BLAST records
    for record in NCBIXML.parse(open(output_file)):
        # Iterate through the alignments
        for alignment in record.alignments:
            # Iterate through the HSPs
            for hsp in alignment.hsps:
                # Check if the HSP meets the e-value threshold
                if hsp.expect > evalue:
                    with open("specificity_results.txt", "a") as f:
                        f.write("query: %s" % record.query[:100])
                        f.write("Match: {}, E-value: {}".format(alignment.title, hsp.expect))
                        f.write("\n")
                        f.write("Coverage: {} %".format((hsp.query_end - hsp.query_start)/len(record.query)*100))
