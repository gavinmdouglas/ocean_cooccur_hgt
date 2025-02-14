#!/usr/bin/python3

import gzip
from collections import defaultdict

# Identify genomes with only one chromosome (in addition to any plasmids).

genome_count = defaultdict(int)

with gzip.open('/mfs/gdouglas/projects/ocean_mags/progenomes_analyses/genome_prep/representatives.aquatic.contigs.headers.txt.gz', 'rt') as fasta_header_fh:
    for line in fasta_header_fh:
        if line.startswith('>'):

            if 'plasmid' in line.lower():
                continue

            line = line[1:].strip()
            contig_subset = line.split(' ')[0]
            contig_subset_split = contig_subset.split('.')
            genome_id = contig_subset_split[0] + '.' + contig_subset_split[1]

            genome_count[genome_id] += 1


with open('/mfs/gdouglas/projects/ocean_mags/progenomes_analyses/genome_prep/single_chrom_genomes.txt', 'w') as outfile:
    for genome_id, chrom_count in genome_count.items():
        if chrom_count == 1:
            outfile.write(genome_id + '\n')
