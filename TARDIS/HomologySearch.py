#!/usr/bin/env python3

# Description
###############################################################################
'''
Protocol implementation to look for specificity in the targets found. 

1. Given a FASTA file with protein sequences, creates a local database
2. Performs a Psi-BLAST search against the database
3. Returns targets with similarity below a given e-value threshold (default: 0.04)

They are imported as:

import TARDIS.HomologySearch as HomologySearch
'''

# Imports
###############################################################################
from TARDIS.Initialise import *

import os
import subprocess

from Bio.Blast import NCBIXML

# Functions
###############################################################################
def make_database(fasta_file: str, database_name: str) -> None:
    '''Make a local database from a fasta file

    Parameters
    ----------
    fasta_file (str)
        Path to the fasta file
    database_name (str)
        Name of the database to be created
    '''

    # make a local database from the fasta file
    make_db_cmd = ["makeblastdb",
                   "-in", fasta_file,
                   "-dbtype", "prot",
                   "-out", database_name,
                   "-title", database_name,
                   "-parse_seqids"]
    
    subprocess.run(make_db_cmd, shell=True)

    return None

def blast_check(database: str, targets: str, evalue: float) -> None:
    '''Perform a Psi-BLAST search against a local database and return targets with similarity below a given e-value threshold.

    Parameters
    ----------
    database (str)
        Path to the database
    targets (str)
        Path to the directory containing the fasta files with the targets
    evalue (float)
        E-value threshold
    '''

    # set the output file name
    output_file = "blast_results.xml"

    # perform psi-blast on all fasta files in the directory
    for out_dir in os.listdir(targets):
        if out_dir.endswith(".fasta"):
            print(f"Blasting {out_dir}")
            query_file = os.path.join(targets, out_dir)
            blast_cmd = [
                "psiblast",
                "-query", query_file,
                "-db", database,
                "-num_iterations", 3,
                "-outfmt", "5"
            ]
            subprocess.run(blast_cmd, shell=True)

    # Iterate through the BLAST records
    for record in NCBIXML.parse(open(output_file)):
        # Iterate through the alignments
        for alignment in record.alignments:
            # Iterate through the HSPs
            for hsp in alignment.hsps:
                # Check if the HSP meets the e-value threshold
                if hsp.expect < evalue:
                    print('\n')
                    print(f"query: {record.query[:100]}")
                    print(f"Match: {alignment.title}, E-value: {hsp.expect}")

    return None
