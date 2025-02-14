#!/usr/bin/env python3
import gzip
from collections import defaultdict

# Read in FASTA file.
fasta = defaultdict(list)
genome_to_seqid = defaultdict(set)

with gzip.open('/mfs/gdouglas/projects/ocean_mags/progenomes_analyses/download/representatives.aquatic.contigs.fasta.gz', 'rt') as fasta_fh:
    for line in fasta_fh:
        if line.startswith('>'):
            sequence_id = line[1:].strip()
            contig_subset = sequence_id.split(' ')[0]
            contig_subset_split = contig_subset.split('.')
            genome_id = 'g' + contig_subset_split[0] + '_' + contig_subset_split[1]

            genome_to_seqid[genome_id].add(sequence_id)
        else:
            fasta[sequence_id].append(line.strip().upper())

# Write out FASTA files.
outfolder = '/mfs/gdouglas/projects/ocean_mags/progenomes_analyses/genome_prep/genome_fastas/'

for genome in genome_to_seqid.keys():
    with open(outfolder + genome + '.fa', 'w') as out_fh:
        for sequence_id in genome_to_seqid[genome]:
            out_fh.write('>' + sequence_id + '\n')
            for line in fasta[sequence_id]:
                out_fh.write(line + '\n')
