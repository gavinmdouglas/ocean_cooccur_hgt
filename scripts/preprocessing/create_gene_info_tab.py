#!/usr/bin/python3

import sys
import gzip
import os

# Create a new table of gene info, including gene representative, scaffold, and genome.
# Ignore all genes on contigs < 5000 bp.

def main():

    # Read in genomes to consider.
    passed_genomes = set()
    with gzip.open('/mfs/gdouglas/projects/ocean_mags/mapfiles/MAGs_to_analyze.txt.gz', 'rt') as genomes_to_consider_fh:
        for genomes_to_consider_line in genomes_to_consider_fh:
            passed_genomes.add(genomes_to_consider_line.rstrip())

    # Read in map of scaffold IDs to genome IDs from gzipped file.
    scaffold_to_genome = {}
    with gzip.open('/mfs/gdouglas/projects/ocean_mags/Sunagawa_dataset/Sunagawa_Archaea_Bacteria_scaffold_to_genome.txt.gz', 'rt') as scaffold_to_genome_fh:
        for scaffold_map_line in scaffold_to_genome_fh:
            scaffold_map_line = scaffold_map_line.rstrip().split()
            scaffold_to_genome[scaffold_map_line[0]] = scaffold_map_line[1]

    # Read in mapfile of gene to gene representative.
    gene_to_rep = {}
    with gzip.open('/mfs/gdouglas/projects/ocean_mags/Sunagawa_dataset/gene-catalog-membership.tsv.gz', 'rt') as gene_to_gene_rep_fh:
        for gene_map_line in gene_to_gene_rep_fh:
            gene_map_line = gene_map_line.rstrip().split()
            gene_to_rep[gene_map_line[0]] = gene_map_line[1]

    # Read through all Prokka GFF files and identify (1) which scaffolds are at least 5000 bp long and (2) which genes are on those scaffolds.
    gff_filepaths = list()
    for file in os.listdir('/mfs/gdouglas/projects/ocean_mags/Sunagawa_dataset/prokka'):
        if file.endswith('.gff'):
            gff_filepaths.append(os.path.join('/mfs/gdouglas/projects/ocean_mags/Sunagawa_dataset/prokka', file))

    failed_genome = 0
    passed_genome = 0
    no_gene_rep = 0
    total_genes = 0

    print('gene' + '\t' + 'gene_rep' + '\t' + 'scaffold' + '\t' + 'genome')
    for gff_filepath in gff_filepaths:
        genome_id = '.'.join(gff_filepath.split('/')[-1].split('.')[:-1])

        if genome_id not in passed_genomes:
            failed_genome += 1
            continue
        else:
            passed_genome += 1

        with open(gff_filepath, 'r') as gff_fh:
            passed_scaffolds = set()
            for line in gff_fh:
                if line.startswith('##sequence-region'):
                    line_split = line.split()
                    if int(line_split[3]) >= 5000:
                        passed_scaffolds.add(line_split[1])
                    continue
                elif line.startswith('#'):
                    continue
                else:
                    line_split = line.rstrip().split('\t')
                    if len(line_split) != 9:
                        continue

                    if line_split[0] in passed_scaffolds:
                        # Parse gene ID.
                        gene_id = line_split[8].split(';')[0].split('=')[1]

                        if line_split[2] != 'CDS':
                            continue

                        total_genes += 1
                        
                        # Skip if no gene rep found, which means it was left out of the clustering for some reason
                        # in the Sunagawa paper (but keep track of how many genes are skipped).
                        if gene_id not in gene_to_rep.keys():
                            no_gene_rep += 1
                            continue

                        # Sanity check that scaffold is in expected genome.
                        if scaffold_to_genome[line_split[0]] != genome_id:
                            sys.exit('Scaffold ' + line_split[0] + ' is not expected to be in this genome ' + genome_id + ', and instead was expected to be in: ' + scaffold_to_genome[line_split[0]])
                        
                        print(gene_id + '\t' + gene_to_rep[gene_id] + '\t' + line_split[0] + '\t' + genome_id)

    print('Number of genomes that failed CheckM completeness/contamination thresholds: ' + str(failed_genome),
            file = sys.stderr)
    print('Number of genomes that passed CheckM completeness/contamination thresholds: ' + str(passed_genome),
            file = sys.stderr)
    print('Number of genes that were skipped because no gene representative was found: ' + str(no_gene_rep) + ' out of ' + str(total_genes) + ' total genes',
            file = sys.stderr)
    
if __name__ == '__main__':
    main()
