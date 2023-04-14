from Bio import SeqIO
from Bio import pairwise2

import os

# Set the directory containing the FASTA files
directory = "/path/to/fasta/files/"

# Set the minimum identity score for matches (in percent)
min_identity_score = 90

# Loop through all files in the directory
for filename in os.listdir(directory):
    # Check if file is a FASTA file
    if filename.endswith(".fasta") or filename.endswith(".fa"):
        # Read in the sequences from the file
        sequences = list(SeqIO.parse(os.path.join(directory, filename), "fasta"))

        # Loop through all pairs of sequences
        for i in range(len(sequences)):
            for j in range(i+1, len(sequences)):
                # Align the sequences
                alignments = pairwise2.align.globalxx(sequences[i].seq, sequences[j].seq)

                # Get the highest scoring alignment
                best_alignment = max(alignments, key=lambda x: x.score)

                # Calculate the identity score
                identity_score = (best_alignment.matches / len(best_alignment.seqA)) * 100

                # Print the filenames and identity score if the score is above the minimum
                if identity_score >= min_identity_score:
                    print(f"{filename}: {sequences[i].id} vs {sequences[j].id}: {identity_score:.2f}% identity")
