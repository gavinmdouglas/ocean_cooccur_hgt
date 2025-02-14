#!/usr/bin/env python3
import gzip
from collections import defaultdict
import sys

def reverse_complement(dna):
   complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
   return ''.join(complement[base] for base in reversed(dna.upper()))

# Run for each genome individually, to speed up process.
focal_genome = sys.argv[1]

if focal_genome is None:
    sys.exit('Provide a genome ID.')
elif not focal_genome.startswith('g') and '_' not in focal_genome:
    sys.exit('Provide a genome ID in the format gXX_YY.')

gene_seq = {}
focal_current = False
with gzip.open('/mfs/gdouglas/projects/ocean_mags/progenomes_analyses/download/representatives.aquatic.genes.fasta.gz', 'rt') as gene_fasta_fh:
    for line in gene_fasta_fh:
        if line.startswith('>'):
            line_split = line[1:].strip().split()
            gene_id = line_split[0]
            gene_id_split = gene_id.split('.')
            genome_id = 'g' + gene_id_split[0] + '_' + gene_id_split[1]
            if genome_id == focal_genome:
                gene_seq[gene_id] = ''
                focal_current = True
            else:
                focal_current = False
        elif focal_current:
            gene_seq[gene_id] += line.strip().upper()
print('Read in genes.', file=sys.stderr)

gene_to_contig = {}
contig_seq = {}
focal_current = False
with gzip.open('/mfs/gdouglas/projects/ocean_mags/progenomes_analyses/download/representatives.aquatic.contigs.fasta.gz', 'rt') as fasta_fh:
    for line in fasta_fh:
        if line.startswith('>'):
            sequence_id = line[1:].strip()
            contig_subset = sequence_id.split(' ')[0]
            contig_subset_split = contig_subset.split('.')
            genome_id = 'g' + contig_subset_split[0] + '_' + contig_subset_split[1]
            if genome_id == focal_genome:
                contig_seq[contig_subset] = ''
                focal_current = True
            else:
                focal_current = False
        elif focal_current:
            contig_seq[contig_subset] += line.strip().upper()
print('Read in contigs.', file=sys.stderr)

# print("gene\tgenome\tgene_length\tcontig")
for gene_id in gene_seq.keys():
    gene_length = len(gene_seq[gene_id])
    contig_matches = set()
    reverse_geneseq = reverse_complement(gene_seq[gene_id])
    for contig_id in contig_seq.keys():
        if gene_seq[gene_id] in contig_seq[contig_id] or reverse_geneseq in contig_seq[contig_id]:
            contig_matches.add(contig_id)
    if len(contig_matches) == 1:
        print(gene_id + '\t' + focal_genome + '\t' + str(gene_length) + '\t' + list(contig_matches)[0])
    elif len(contig_matches) > 1:
        contig_matches = sorted(list(contig_matches))
        print(gene_id + '\t' + focal_genome + '\t' + str(gene_length) + '\t' + ';'.join(contig_matches))
    else:
        print(gene_id + '\t' + focal_genome + '\t' + str(gene_length) + '\tNot_found')
