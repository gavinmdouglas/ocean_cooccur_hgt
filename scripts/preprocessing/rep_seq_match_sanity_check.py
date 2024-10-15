#!/usr/bin/python3

import os
import sys

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from functions import read_fasta

# Sanity check that representative CDS sequences (from Sunagawa resource)
# are all present and are exact matches in the newly parsed FASTAs.

rep_seqs = read_fasta('/mfs/gdouglas/projects/ocean_mags/Sunagawa_dataset/gene-catalog.fasta.gz')

parsed_seqs = read_fasta('/mfs/gdouglas/projects/ocean_mags/Sunagawa_dataset/all_cds_sequences.fasta')

for rep_seq_id in rep_seqs:
    if rep_seq_id not in parsed_seqs:
        print(f"Representative sequence {rep_seq_id} not found in parsed sequences.")
    elif rep_seqs[rep_seq_id] != parsed_seqs[rep_seq_id]:
        print(f"Representative sequence {rep_seq_id} does not match parsed sequence.")

print("Sanity check complete.")
