#!/usr/bin/python3

import os
import sys
import gzip
from collections import defaultdict

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from functions import read_fasta

# Confirmation that some genes with the same representative sequences are of quite different lengths.

parsed_seqs = read_fasta('/mfs/gdouglas/projects/ocean_mags/Sunagawa_dataset/all_cds_sequences.fasta')

rep_gene_lengths = defaultdict(list)

missing_gene = 0

with gzip.open('/mfs/gdouglas/projects/ocean_mags/Sunagawa_dataset/gene-catalog-membership.tsv.gz', 'rt') as f:
    header = f.readline().strip().split('\t')
    gene_to_rep_seq = {}
    for line in f:
        parts = line.strip().split('\t')
        gene_id = parts[0]
        rep_id = parts[1]

        if gene_id in parsed_seqs:
            rep_gene_lengths[rep_id].append(len(parsed_seqs[gene_id]))
        else:
            missing_gene += 1

num_major_diff = 0
no_major_diff = 0
singleton = 0
for rep_id, lengths in rep_gene_lengths.items():
    if len(lengths) > 1:
        # Check if the min and max lengths differ by more than 50%.
        min_len = min(lengths)
        max_len = max(lengths)
        if max_len > 1.5 * min_len:
            num_major_diff += 1
            print(f'Representative ID: {rep_id}, min: {min_len}, max: {max_len}')
        else:
            no_major_diff += 1
    else:
        singleton += 1

print('\n\n\n')
print(f'Number of representative sequences with major length differences: {num_major_diff}')
print(f'Number of representative sequences without major length differences: {no_major_diff}')
print(f'Number of singleton representative sequences: {singleton}')
print(f'Number of genes missing from parsed sequences: {missing_gene}')