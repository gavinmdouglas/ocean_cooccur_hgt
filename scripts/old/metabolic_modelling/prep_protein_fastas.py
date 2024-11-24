#!/usr/bin/python3

import gzip
import os
import sys
from collections import defaultdict

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from functions import write_fasta, read_fasta

# Parse all protein FASTAs and get proteins per genome.
# Note that this is only for the genomes that are metabolically modelled,
# which corresponds to the set in the allsamples presence file.

# Read in allsamples presence file to get taxa of interest.
presence_file = "/mfs/gdouglas/projects/ocean_mags/coverm/combined_tables/metaG_presence_allsamples.tsv.gz"
with gzip.open(presence_file, 'rt') as f:
    taxa = set(f.readline().strip().split("\t")[1:])

# Read in mapping from taxon to genome.
with gzip.open("/mfs/gdouglas/projects/ocean_mags/mapfiles/MAG_to_taxon.tsv.gz", 'rt') as taxon_map_fh:
    taxon_map = dict()
    for line in taxon_map_fh:
        line_split = line.split()
        if line_split[1] in taxa:
            taxon_map[line_split[0]] = line_split[1]

genomes = set(list(taxon_map.keys()))

# Read mapping of gene ID to genome.
gene_map = {}
with open("/mfs/gdouglas/projects/ocean_mags/Sunagawa_dataset/gene_info_allscaffolds.tsv", 'r') as gene_map_fh:
    gene_map_fh.readline()
    for gene_line in gene_map_fh:
        gene_line = gene_line.strip().split()
        gene_map[gene_line[0]] = gene_line[3]

# Read map of scaffolds to genomes.
scaffold_to_genome = {}
with gzip.open("/mfs/gdouglas/projects/ocean_mags/Sunagawa_dataset/Sunagawa_Archaea_Bacteria_scaffold_to_genome.txt.gz", 'rt') as scaffold_to_genome_fh:
    for line in scaffold_to_genome_fh:
        line_split = line.strip().split()
        scaffold_to_genome[line_split[0]] = line_split[1]

missing_scaffolds = set()

# Read in all protein sequences and allocate to appropriate genomes.
all_proteins = read_fasta("/mfs/gdouglas/projects/ocean_mags/Sunagawa_dataset/gene-catalog.faa")
genome_seqs = defaultdict(dict)
present_genes = 0
present_focal_genes = 0
for name_raw, seq in all_proteins.items():
    # Remove trailing "_1" from name.
    name = name_raw[:-2]
    name_split = name.split("-")
    scaffold = "-".join(name_split[:-1])

    if scaffold not in scaffold_to_genome.keys():
        missing_scaffolds.add(scaffold)
        continue

    genome = scaffold_to_genome[scaffold]
    if genome not in genomes:
        continue

    if name not in gene_map.keys():
        sys.exit("Gene ID " + name + " not found in gene map.")
    else:
        present_genes += 1

    if genome != gene_map[name]:
        sys.exit("Error: Gene " + name + " not found in genome " + genome + ".")

    if genome not in genomes:
        continue
    elif name in genome_seqs[genome].keys():
        sys.exit("Error: protein ID " + name + " already found in genome " + genome + ".")
    else:
        genome_seqs[genome][name] = seq
        present_focal_genes += 1

print("Missing scaffolds: " + str(len(missing_scaffolds)))
print("Present genes: " + str(present_genes))
print("Present focal genes: " + str(present_focal_genes))

for genome in genome_seqs:
    taxon = taxon_map[genome]
    num_proteins = len(genome_seqs[genome])
    write_fasta(genome_seqs[genome], "/mfs/gdouglas/projects/ocean_mags/metabolic_modelling/protein_fastas/" + taxon + ".faa")
