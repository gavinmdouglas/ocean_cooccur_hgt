#!/usr/bin/python3

import gzip
import os
import sys
from collections import defaultdict

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from functions import write_fasta, read_fasta

gene_map = set()
with gzip.open("/mfs/gdouglas/projects/ocean_mags/Sunagawa_dataset/gene_info.tsv.gz", 'rt') as gene_map_fh:
    gene_map_fh.readline()
    for gene_line in gene_map_fh:
        gene_line = gene_line.strip().split()
        if gene_line[3] == "TARA_SAMEA2622518_METAG_FNCNCLIF":
            gene_map.add(gene_line[0])

all_proteins = read_fasta("/mfs/gdouglas/projects/ocean_mags/Sunagawa_dataset/gene-catalog.faa")

present = 0
missing = 0
for gene_name in gene_map.keys():
    if gene_name in all_proteins:
        present += 1
    else:
        missing += 1

print("Present: " + str(present))
print("Missing: " + str(missing))
