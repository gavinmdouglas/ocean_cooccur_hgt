#!/usr/bin/python3

import gzip
from collections import defaultdict

# Get tally of separate contigs, by whether they are contigs or plasmids.

plasmid_count = defaultdict(int)
other_count = defaultdict(int)

with gzip.open('/mfs/gdouglas/projects/ocean_mags/progenomes_analyses/genome_prep/representatives.aquatic.contigs.headers.txt.gz', 'rt') as fasta_header_fh:
    for line in fasta_header_fh:
        if line.startswith('>'):
            line = line[1:].strip()
            contig_subset = line.split(' ')[0]
            contig_subset_split = contig_subset.split('.')
            genome_id = contig_subset_split[0] + '.' + contig_subset_split[1]

            if 'plasmid' in line.lower():
                plasmid_count[genome_id] += 1
            else:
                other_count[genome_id] += 1

unique_genomes = sorted(set(list(plasmid_count.keys()) + list(other_count.keys())))
with open('/mfs/gdouglas/projects/ocean_mags/progenomes_analyses/genome_prep/contig_tally_plasmid_vs_other.tsv', 'w') as outfile:
    outfile.write('genome_id\tplasmid_count\tother_count\n')
    for genome_id in unique_genomes:
        outfile.write(f'{genome_id}\t{plasmid_count[genome_id]}\t{other_count[genome_id]}\n')
