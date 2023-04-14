import os
import subprocess
from Bio.Blast import NCBIXML

# set the path to the local blast database
db = '/home/camila/LMDM/Mestrado/TARDIS/TARDIS/humanProtDB/humanProtDB'

# set the directory containing the fasta files
targets = '/home/camila/LMDM/Mestrado/TARDIS/ESKAPE/targets_Pseudomonas_aeruginosa'

# set the e-value threshold
evalue_thresh = 0.04

# set the output file name
output_file = "blast_results.xml"

# perform blastp on all fasta files in the directory
for out_dir in os.listdir(targets):
    if out_dir.endswith('.fasta'):
        print("Blasting " + out_dir)
        query_file = os.path.join(targets, out_dir)
        cmd = "blastp -query " + query_file + " -db " + db + " -out " + output_file + " -evalue 0.001 -outfmt 5"
        subprocess.run(cmd, shell=True)


# Iterate through the BLAST records
for record in NCBIXML.parse(open(output_file)):
    # Iterate through the alignments
    for alignment in record.alignments:
        # Iterate through the HSPs
        for hsp in alignment.hsps:
            # Check if the HSP meets the e-value threshold
            if hsp.expect < evalue_thresh:
                print("\n")
                print("query: %s" % record.query[:100])
                print("Match: {}, E-value: {}".format(alignment.title, hsp.expect))
