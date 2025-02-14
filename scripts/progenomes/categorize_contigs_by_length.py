#!/usr/bin/env python3
import gzip
import sys
from collections import defaultdict

# Read in FASTA file.
fasta = defaultdict(str)
genome_to_seqid = defaultdict(set)
seqid_to_genome = {}

plasmid_seqids = set()

length_cutoff = 10000

with gzip.open('/mfs/gdouglas/projects/ocean_mags/progenomes_analyses/download/representatives.aquatic.contigs.fasta.gz', 'rt') as fasta_fh:
    for line in fasta_fh:
        if line.startswith('>'):
            sequence_id = line[1:].strip()
            contig_subset = sequence_id.split(' ')[0]
            contig_subset_split = contig_subset.split('.')
            genome_id = 'g' + contig_subset_split[0] + '_' + contig_subset_split[1]
            genome_to_seqid[genome_id].add(contig_subset)
            seqid_to_genome[contig_subset] = genome_id
            if 'plasmid' in sequence_id or 'Plasmid' in sequence_id:
                plasmid_seqids.add(contig_subset)
        else:
            fasta[contig_subset] += line.strip()

print("Done reading in FASTA.", file=sys.stderr)

genome_long_counts = defaultdict(int)
genome_short_counts = defaultdict(int)

with open('/mfs/gdouglas/projects/ocean_mags/progenomes_analyses/genome_prep/contig_length_breakdown.tsv', 'w') as contig_breakdown_fh:
    contig_breakdown_fh.write("genome\tcontig\ttype\tlength\tplasmid\n")
    for genome_id, seqids in genome_to_seqid.items():
        for seqid in seqids:
            sequence = fasta[seqid]
            seq_len = len(sequence)

            if seqid in plasmid_seqids:
                plasmid_info = 'Plasmid'
            else:
                plasmid_info = 'Not_plasmid_labelled'

            if seq_len >= length_cutoff:
                genome_long_counts[genome_id] += 1
                contig_breakdown_fh.write(genome_id + '\t' + seqid + '\t' + 'long' + '\t' + str(seq_len) + '\t' + plasmid_info + '\n')
            else:
                genome_short_counts[genome_id] += 1
                contig_breakdown_fh.write(genome_id + '\t' + seqid + '\t' + 'short' + '\t' + str(seq_len) + '\t' + plasmid_info + '\n')

with open('/mfs/gdouglas/projects/ocean_mags/progenomes_analyses/genome_prep/genome_length_breakdown.tsv', 'w') as genome_breakdown_fh:
    genome_breakdown_fh.write('genome\tlong_contigs\tshort_contigs\n')
    for genome_id in sorted(list(genome_to_seqid.keys())):
        genome_breakdown_fh.write(genome_id + '\t' + str(genome_long_counts[genome_id]) + '\t' + str(genome_short_counts[genome_id]) + '\n')
