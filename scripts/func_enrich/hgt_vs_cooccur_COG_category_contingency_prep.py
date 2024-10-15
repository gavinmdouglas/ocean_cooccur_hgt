#!/usr/bin/python3

import argparse
import pandas as pd
from collections import defaultdict
import gzip
import sys


def main():

    parser = argparse.ArgumentParser(

        description='''
Prep contingency tables for COG categories enriched among BLAST hits between co-occurring and non-cooccuring genomes.
By taxonomic level and identity cut-off separately.
''',

        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('--cooccur_tab', metavar='COOCCUR_TAB', type=str,
                        help="Path to gzipped co-occurrence table. Must be tab. delimited and start with the columns taxon_i and taxon_j",
                        required=True)

    parser.add_argument('--coverm_tab', metavar="MEASURE", type=str,
                        help="Column(s) to use from co-occurrence table to use as association measure. Input comma-delimited commas if multiple desired.",
                        required=True)

    parser.add_argument('--cooccur_measure', metavar="MEASURE", type=str,
                        help="Column name to use for co-occurrence measure.",
                        required=False,
                        default="ratio")

    parser.add_argument('--cooccur_measure_p', metavar="MEASURE", type=str,
                        help="Column name to use for co-occurrence measure P-value.",
                        required=False,
                        default="BH")

    parser.add_argument('--gene_beds_folder', metavar='PATH', type=str,
                        help="Path to folder containing gene hit bedfiles.",
                        required=False,
                        default= "/mfs/gdouglas/projects/ocean_mags/blast_output/intersections/blast_hit_gene_beds/")

    parser.add_argument('--cooccur_cutoff', metavar="MEASURE", type=float,
                        help="Estimate cut-off to use to identify a significant hit as higher or lower than expected by chance.",
                        required=False,
                        default=1.0)

    parser.add_argument("--taxa_tab", metavar="TAX_TAB", type=str,
                        help="Path to taxonomic table breakdown (gzipped).",
                        required=False,
                        default="/mfs/gdouglas/projects/ocean_mags/mapfiles/MAG_taxa_breakdown.tsv.gz")

    parser.add_argument("--passing_scaffolds", metavar="PATH", type=str,
                        help="Path to file with all passing scaffolds to consider.",
                        required=False,
                        default="/mfs/gdouglas/projects/ocean_mags/Sunagawa_dataset/scaffolds_5000bp.txt")

    parser.add_argument("--full_gene_info", metavar="PATH", type=str,
                        help="Path to (gzipped) table with information on each gene (gene rep., scaffold, and genome).",
                        required=False,
                        default="/mfs/gdouglas/projects/ocean_mags/Sunagawa_dataset/gene_info.tsv.gz")

    parser.add_argument("--COG_annot", metavar="PATH", type=str,
                        help="Path to (gzipped) table of COG annotations per representative gene.",
                        required=False,
                        default="/mfs/gdouglas/projects/ocean_mags/Sunagawa_dataset/genomes-representative-cog-info.tsv.gz")

    args = parser.parse_args()

    higher_levels = ["Domain", "Phylum", "Class", "Order", "Family", "Genus"]

    # Figure out taxa to consider (i.e., those in CoverM table).
    with gzip.open(args.coverm_tab, 'rt') as coverm_in:
        header = coverm_in.readline().strip().split('\t')
        if header[0] != 'sample':
            sys.exit("First column of coverm table must be 'sample'")
        focal_taxa = set(header[1:])

    # Read in mapping of genome ID to taxonomic category for each taxonomic levels.
    MAG_to_taxon = dict()
    taxa_map = defaultdict(dict)
    taxa_tab = pd.read_table(filepath_or_buffer=args.taxa_tab, sep='\t',
                             compression='gzip', header=0, index_col=1, na_values=['NA'])
    for taxon in taxa_tab.index.values:
        if taxon not in focal_taxa:
            continue

        for level in higher_levels:
            taxa_map[taxon][level] = taxa_tab.loc[taxon, level]

        MAG_to_taxon[taxa_tab.loc[taxon, 'MAG']] = taxon

    # Identify sig. co-occurring MAG pairs.
    cooccur_tab = pd.read_table(args.cooccur_tab, compression='gzip', header=0, index_col=None)
    positive_cooccur_combos = set()
    for index, row in cooccur_tab.iterrows():
        if row[args.cooccur_measure] > args.cooccur_cutoff and row[args.cooccur_measure_p] < 0.05:
            positive_cooccur_combos.add(",".join([row["taxon_i"], row["taxon_j"]]))
            positive_cooccur_combos.add(",".join([row["taxon_j"], row["taxon_i"]]))

    # Read in passing scaffolds from file.
    passing_scaffolds = set()
    with open(args.passing_scaffolds, 'r') as f:
        for line in f:
            passing_scaffolds.add(line.strip())

    # Read in full gene information.
    gene_to_genome = dict()
    gene_to_rep = dict()
    scaffold_to_genome = dict()
    with gzip.open(args.full_gene_info, 'rt') as info_fh:
        info_fh.readline()
        for info_line in info_fh:
            gene, gene_rep, scaffold, genome = info_line.strip().split('\t')
            if scaffold not in passing_scaffolds:
                continue
            gene_to_genome[gene] = genome
            gene_to_rep[gene] = gene_rep
            scaffold_to_genome[scaffold] = genome

    # Read in COG annotations.
    all_categories = set()
    rep_to_COG_categories = dict()
    with gzip.open(args.COG_annot, 'rt') as COG_fh:
        COG_fh.readline()
        for COG_line in COG_fh:
            rep_gene, COG, COG_category_raw = COG_line.strip().split('\t')

            if COG_category_raw == 'NA' or COG_category_raw == '' or COG_category_raw == '-':
                continue

            rep_to_COG_categories[rep_gene] = COG_category_raw.split(',')
            all_categories.update(set(rep_to_COG_categories[rep_gene]))

    all_categories = sorted(list(all_categories))

    print("Read in initial input files", file=sys.stderr)

    cutoffs = ['95', '99']

    # Print header.
    print("Identity\tTax_Level\tCOG_Category\tCooccur_HGT_Gene_Count\tNoncooccur_HGT_Gene_Count")

    gene_missing = 0
    scaffold2_missing = 0

    for cutoff in cutoffs:

        print(f"Processing for identity cut-off: {cutoff}", file=sys.stderr)

        for tax_level in higher_levels:

            print(f"Taxa level {tax_level}", file=sys.stderr)

            cooccur_COG_breakdown = defaultdict(int)
            noncooccur_COG_breakdown = defaultdict(int)

            gene_bed_file = args.gene_beds_folder + f"/Hits_{tax_level}_{cutoff}_genes_with_B.bed"
            with open(gene_bed_file, 'r') as gene_bed_fh:
                for line in gene_bed_fh:
                    line_split = line.strip().split('\t')
                    scaffold1 = line_split[0]
                    scaffold2 = line_split[13].split('||')[0].replace('subject_', '')

                    if scaffold1 not in passing_scaffolds or scaffold2 not in passing_scaffolds:
                        continue

                    gene = line_split[3]

                    if gene not in gene_to_genome:
                        gene_missing += 1
                        continue

                    if scaffold2 not in scaffold_to_genome:
                        scaffold2_missing += 1
                        continue

                    genome1 = gene_to_genome[gene]
                    genome2 = scaffold_to_genome[scaffold2]
                    if genome1 not in MAG_to_taxon or genome2 not in MAG_to_taxon:
                        continue

                    gene_rep = gene_to_rep[gene]
                    if gene_rep not in rep_to_COG_categories:
                        continue

                    taxon1 = MAG_to_taxon[genome1]
                    taxon2 = MAG_to_taxon[genome2]

                    cooccuring_flag = False
                    if ",".join([taxon1, taxon2]) in positive_cooccur_combos:
                        cooccuring_flag = True
                    elif ",".join([taxon2, taxon1]) in positive_cooccur_combos:
                        cooccuring_flag = True

                    if cooccuring_flag:
                        for COG_category in rep_to_COG_categories[gene_rep]:
                            cooccur_COG_breakdown[COG_category] += 1
                    else:
                        for COG_category in rep_to_COG_categories[gene_rep]:
                            noncooccur_COG_breakdown[COG_category] += 1

            for category in all_categories:
                print(f"{cutoff}\t{tax_level}\t{category}\t{cooccur_COG_breakdown[category]}\t{noncooccur_COG_breakdown[category]}")

    print(f"Number of genes missing from gene info table: {gene_missing}", file=sys.stderr)
    print(f"Number of scaffolds missing from scaffold to genome table: {scaffold2_missing}", file=sys.stderr)

if __name__ == '__main__':
    main()
