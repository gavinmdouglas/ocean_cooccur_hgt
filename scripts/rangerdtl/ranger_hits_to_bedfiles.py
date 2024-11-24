#!/usr/bin/python3

import argparse
from collections import defaultdict
import gzip
import sys


def main():

    parser = argparse.ArgumentParser(

            description='''
            Parse simplified RANGER-DTL hits breakdown, along with gene annotation information, and write out
            two bedfiles: (1) all putatively horizontally transferred genes separately and (2) all other genes.
            Will ignore genes on scaffolds below a certain length threshold (default: 5000 bp).
            Assumes that GFFs are named by genome ID indicated in the simplified RANGER-DTL file and end in ".gff".
            Note that Panaroo-identified 'refound' genes are ignored.
            Can optionally include all RANGER simplified hits (default), or direct hits only (with --direct_only).
            ''',

            epilog='''Usage example:

    python ranger_hits_to_bedfiles.py --ranger TAB --prokka_gff GFF --panaroo TAB --direct_only --output OUTPUT

    ''', formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-r', '--ranger', metavar='INPUT', type=str, required=True,
                        help='Simplified RANGER-DTL summary of all putative HGT events (e.g., simplified_pairwise_transfers.tsv).')

    parser.add_argument('--prokka_gff_folder', metavar='FILE', type=str, required=True,
                        help='Path to folder containing Prokka GFFs.')

    parser.add_argument('--panaroo', metavar='FILE', type=str, required=True,
                        help='Path to Panaroo pangenome presence/absence table to parse.')

    parser.add_argument('--genomes', metavar='FILE', type=str, required=True,
                        help='Path files containing genome IDs to consider (transfers between any other genomes will be ignored). This can include genomes from other species that are not included here.')

    parser.add_argument('--direct_only', action='store_true',
                        help='Only include direct hits (as defined in the simplified RANGER output) as HGT-associated genes.',
                        required=False)

    parser.add_argument('-m', '--min_length', metavar='INT', type=int, required=False,
                        default=5000,
                        help='Path to Panaroo pangenome presence/absence table to parse.')

    parser.add_argument('-o', '--output_prefix', metavar='OUTPUT', type=str, required=True,
                        help='Prefix for two output files (which will be in BED format).')

    args = parser.parse_args()

    # Read through file with genome IDs to consider.
    genome_ids_to_consider = set()
    with open(args.genomes, 'r') as genomes_filehandle:
        for line in genomes_filehandle:
            genome_ids_to_consider.add(line.rstrip())

    # Read through RANGER-DTL file once to get all unique genome IDs.
    genome_ids = set()
    with open(args.ranger, 'r') as ranger_filehandle:
        ranger_filehandle.readline()
        for line in ranger_filehandle:
            line_split = line.rstrip().split('\t')
            genome1 = line_split[1]
            genome2 = line_split[2]

            if genome1 in genome_ids_to_consider and genome2 in genome_ids_to_consider:
                genome_ids.add(genome1)
                genome_ids.add(genome2)

    # Read through each GFF file and get mapping of gene IDs to scaffold IDs.
    # Also keep track of scaffolds above length cut-off.
    # Last, for all genes in passing scaffolds, get prepped BED format coor.
    gene_to_scaffold = dict()
    gene_to_genome = dict()
    passed_scaffolds = set()
    genes_bed_format = dict()

    # Identify all filenames in specified folder with extension ".gff".
    for genome_id in sorted(list(genome_ids)):
        gff_filename = args.prokka_gff_folder + '/' + genome_id + '.gff'
        with open(gff_filename, 'r') as gff_filehandle:
            for gff_line in gff_filehandle:

                if gff_line.startswith('##sequence-region'):
                    gff_line_split = gff_line.split()
                    if int(gff_line_split[3]) >= args.min_length:
                        passed_scaffolds.add(gff_line_split[1])
                    continue

                gff_line_split = gff_line.split('\t')

                if len(gff_line_split) < 9 or gff_line_split[2] != 'CDS':
                    continue

                scaffold = gff_line_split[0]
                gene = gff_line_split[8].split(';')[0].replace('ID=', '')

                gene_to_scaffold[gene] = scaffold

                if scaffold in passed_scaffolds:
                    genes_bed_format[gene] = [scaffold, int(gff_line_split[3]) - 1, int(gff_line_split[4]), gene, '.', gff_line_split[6]]
                    gene_to_genome[gene] = genome_id

        # Read through Panaroo presence/absence table and get mapping of gene IDs to gene family IDs.
        genefamily_to_gene = defaultdict(dict)
        gene_to_genefamily = dict()
        with gzip.open(args.panaroo, 'rt') as panaroo_filehandle:
            panaroo_header = panaroo_filehandle.readline()
            panaroo_genomes = panaroo_header.rstrip().split(',')[3:]

            for panaroo_line in panaroo_filehandle:
                panaroo_line_split = panaroo_line.rstrip().split(',')
                gene_family = panaroo_line_split[0]
                genes = panaroo_line_split[3:]
                for i, gene in enumerate(genes):
                    if gene == '':
                        continue
                    genefamily_to_gene[panaroo_genomes[i]][gene_family] = gene

                    # Also get mapping of individual genes to gene families.
                    for g in gene.split(';'):
                        if 'refound' in g:
                            continue

                        if g.endswith('_len'):
                            g = g[:-4]
                        
                        gene_to_genefamily[g] = gene_family

    # Then read through simplified RANGER-DTL output and identify genes to write to each BED file.
    hgt_associated_genes = set()

    with open(args.ranger, 'r') as ranger_filehandle:
        ranger_filehandle.readline()
        for line in ranger_filehandle:
            line_split = line.rstrip().split('\t')
            gene_family = line_split[0]
            genome1 = line_split[1]
            genome2 = line_split[2]

            # Skip if both genomes are not in the list of genomes to consider.
            if genome1 not in genome_ids_to_consider or genome2 not in genome_ids_to_consider:
                continue

            # Skip if only directly inferred transfers should be considered.
            if args.direct_only and line_split[3] != 'direct':
                continue

            gene1_split = genefamily_to_gene[genome1][gene_family].split(';')
            gene2_split = genefamily_to_gene[genome2][gene_family].split(';')

            for gene1 in gene1_split:
                if 'refound' in gene1:
                    continue

                if gene1.endswith('_len'):
                    gene1 = gene1[:-4]

                scaffold1 = gene_to_scaffold[gene1]

                for gene2 in gene2_split:
                    if 'refound' in gene2:
                        continue

                    if gene2.endswith('_len'):
                        gene2 = gene2[:-4]

                    scaffold2 = gene_to_scaffold[gene2]

                    if scaffold1 in passed_scaffolds and scaffold2 in passed_scaffolds:
                        hgt_associated_genes.add(gene1)
                        hgt_associated_genes.add(gene2)

    # Then loop through all genes and write them to the appropriate
    # BED file (HGT-associated or not).

    output_hgt_fh = open(args.output_prefix + '.hgt.bed', 'w')
    output_nonhgt_fh = open(args.output_prefix + '.nonhgt.bed', 'w')
    panaroo_missing_gene = 0
    for gene in sorted(list(genes_bed_format.keys())):

        outline_raw = genes_bed_format[gene]

        if gene not in gene_to_genefamily:
            panaroo_missing_gene += 1
            continue
        outline_raw[3] = gene_to_genome[gene] + '||' + gene_to_genefamily[gene] + '||' + gene
        outline_raw = [str(x) for x in outline_raw]
        outline = '\t'.join(outline_raw)

        if gene in hgt_associated_genes:
            print(outline, file=output_hgt_fh)
        else:
            print(outline, file=output_nonhgt_fh)

    output_hgt_fh.close()
    output_nonhgt_fh.close()

    print('Number of genes in GFFs, but not in Panaroo presence/absence table: ' +
          str(panaroo_missing_gene) +
          ' (out of ' +
          str(len(genes_bed_format)) +
          ')',
          file=sys.stderr)

if __name__ == '__main__':
    main()
