import os
import csv
import parasail

# Set the directory containing the FASTA files
Acineto_dir = "/home/camila/LMDM/TARDIS/ESKAPE/targets_Acinetobacter_baumannii"
Pseudomonas_dir = "/home/camila/LMDM/TARDIS/ESKAPE/targets_Pseudomonas_aeruginosa"
Staphylococcus_dir = "/home/camila/LMDM/TARDIS/ESKAPE/targets_Staphylococcus_aureus"
Streptococcus_dir = "/home/camila/LMDM/TARDIS/ESKAPE/targets_Streptococcus_pneumoniae"
Enterococcus_dir = "/home/camila/LMDM/TARDIS/ESKAPE/targets_Enterococcus_faecalis"
work_dir = "/home/camila/LMDM/Mestrado/TARDIS/ESKAPE"

# Get a list of all the FASTA files in the directory
fasta_files = [f for f in os.listdir(fasta_dir) if f.endswith(".fasta")]

#Dictionary with sequences and their respective directories
seqs = {}
for f in os.listdir(work_dir):
    if os.path.isdir(os.path.join(work_dir,f)):
        dir = os.path.join(work_dir, f)
    for file in os.listdir(dir):
        if file.endswith(".fasta"):
            with open(os.path.join(dir, file), "r") as seq_file:
                seq = seq_file.read()
                seqs[seq] = [str(dir)]
                
                #seqs[str(seq_file.readlines())] = str(dir)
                
     
print(seqs.items())

# Set the alignment parameters
matrix = parasail.blosum62
gap_open = 10
gap_extend = 1

# Initialize a list to store the results
results = []

# Loop over all pairs of FASTA files and perform pairwise alignments
for i in range(len(fasta_files)):
    for j in range(i+1, len(fasta_files)):
        # Read in the sequences from the FASTA files
        with open(os.path.join(fasta_dir, fasta_files[i]), "r") as f1:
            seq1 = f1.read().strip()
        with open(os.path.join(fasta_dir, fasta_files[j]), "r") as f2:
            seq2 = f2.read().strip()
        # Perform the pairwise alignment
        alignment = parasail.sg_stats_striped_16(seq1, seq2, gap_open, gap_extend, matrix)
        # Calculate the identity score in percentage
        percent_id = (alignment.matches)/alignment.length*100
        # Add the results to the list
        results.append((fasta_files[i], fasta_files[j], percent_id))

# Sort the results by identity score from highest to lowest
results = sorted(results, key=lambda x: x[2], reverse=True)

# Write the results to a CSV file
with open("alignment_results.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["File 1", "File 2", "Identity Score"])
    for r in results:
        writer.writerow([r[0], r[1], r[2]])
        print(f"{r[0]}\t{r[1]}\t{r[2]}")
