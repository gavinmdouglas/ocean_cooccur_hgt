#!/usr/bin/python3

import gzip
import argparse
import sys
import pandas as pd
from itertools import combinations


def main():

    parser = argparse.ArgumentParser(

        description='''
Parse co-occurrence table, tips distances, MAG taxa levels, and HGT summary.
Output a table (to STDOUT) with the key information from these tables for each tip pair present.
''',

        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('--cooccur_measure', metavar="MEASURE", type=str,
                        help="Path to gzipped, tab-delimited coverm table (either presence or RPKM). Used to get the set of taxa to consider. I.e., only the header-line is parsed.",
                        required=True)

    parser.add_argument('--cooccur_tab', metavar='COOCCUR_TAB', type=str,
                        help="Path to gzipped co-occurrence table. Must be tab. delimited and start with the columns taxon_i and taxon_j",
                        required=True)

    parser.add_argument('--coverm_tab', metavar="MEASURE", type=str,
                        help="Column(s) to use from co-occurrence table to use as association measure. Input comma-delimited commas if multiple desired.",
                        required=True)

    parser.add_argument("--tip_distances", metavar="TIP_DIST", type=str,
                        help="Path to tip distances table (gzipped).",
                        required=False,
                        default="/mfs/gdouglas/projects/ocean_mags/phylogenetic_analyses/tip_dist.tsv.gz")

    parser.add_argument("--hgt_tab", metavar="HGT_TAB", type=str,
                        help="Path to HGT pairwise tallies (gzipped).",
                        required=False,
                        default="/mfs/gdouglas/projects/ocean_mags/blast_output/pairwise_tally_summary.tsv.gz")

    parser.add_argument("--taxa_tab", metavar="TAX_TAB", type=str,
                        help="Path to taxonomic table breakdown (gzipped).",
                        required=False,
                        default="/mfs/gdouglas/projects/ocean_mags/mapfiles/MAG_taxa_breakdown.tsv.gz")

    args = parser.parse_args()

    # Read in required tables.
    cooccur_tab = pd.read_table(filepath_or_buffer=args.cooccur_tab, compression='gzip', header=0, index_col=None)
    # Get index for co-occur tab, based on sorted taxon_i and taxon_j values per row.

    # If first two column names of cooccur_tab are taxoni and taxonj, then change to taxon_i, and taxon_j.
    if cooccur_tab.columns[0] == "taxoni" and cooccur_tab.columns[1] == "taxonj":
        cooccur_tab.columns = ["taxon_i", "taxon_j"] + list(cooccur_tab.columns[2:])

    cooccur_tab.index = cooccur_tab.apply(lambda row: ",".join(sorted([row["taxon_i"], row["taxon_j"]])), axis=1)
    tip_dist = pd.read_table(filepath_or_buffer=args.tip_distances, compression='gzip', header=0, index_col=0)
    taxa_tab = pd.read_table(filepath_or_buffer=args.taxa_tab, compression='gzip', header=0, index_col=1)
    hgt_tab = pd.read_table(filepath_or_buffer=args.hgt_tab, compression='gzip', header=0, index_col=0)

    with gzip.open(args.coverm_tab, 'rt') as coverm_in:
        header = coverm_in.readline().strip().split('\t')
        if header[0] != 'sample':
            sys.exit("First column of coverm table must be 'sample'")
        taxa = sorted(header[1:])
    taxa_combinations = list(combinations(taxa, 2))

    measures = args.cooccur_measure.split(',')
    if len(measures) == 0:
        sys.exit("No co-occurrence measures provided.")
    headerline = ["taxa_combo", "taxon_i", "taxon_j", "tip_dist", "diff_tax_level", "95_hit_count", "99_hit_count", "both_hit_count", "95_gene_count", "99_gene_count", "both_gene_count"]
    for measure in measures:
        headerline.append("cooccur_" + measure)
    print("\t".join(headerline))

    for taxa_pair in taxa_combinations:
        taxon_i = taxa_pair[0]
        taxon_j = taxa_pair[1]
        taxa_combo = taxon_i + ',' + taxon_j
        reverse_taxa_combo = taxon_j + ',' + taxon_i

        # Check that both taxa IDs are in the tip distance index.
        if taxon_i not in tip_dist.index or taxon_j not in tip_dist.index:
            tip_dist_val = "NA"
        else:
            tip_dist_val = tip_dist.loc[taxon_i, taxon_j]

        diff_level = "NO_DIFF_FOUND"
        tax_levels = ["Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain"]
        taxon_i_levels = taxa_tab.loc[taxon_i, tax_levels]
        taxon_j_levels = taxa_tab.loc[taxon_j, tax_levels]
        for i in range(len(taxon_i_levels)):
            if taxon_i_levels.iloc[i] != taxon_j_levels.iloc[i]:
                diff_level = tax_levels[i]
                break

        # If taxa_combo is in hgt table, then specify values for "95_hit_count", "99_hit_count", "both_hit_count", "95_gene_count", "99_gene_count", "both_gene_count".
        if taxa_combo in hgt_tab.index:
            hgt_row = hgt_tab.loc[taxa_combo]
            hgt_values = [str(hgt_row["95_hit_count"]), str(hgt_row["99_hit_count"]), str(hgt_row["both_hit_count"]), str(hgt_row["95_gene_count"]), str(hgt_row["99_gene_count"]), str(hgt_row["both_gene_count"])]
        elif reverse_taxa_combo in hgt_tab.index:
            sys.exit(f"Error - Reverse taxa_combo found in HGT table: {reverse_taxa_combo}.")
        else:
            hgt_values = ["0"] * 6

        # Get co-occurrence values for each measure.
        cooccur_values = []
        if taxa_combo in cooccur_tab.index:
            cooccur_row = cooccur_tab.loc[taxa_combo]
            for measure in measures:
                cooccur_values.append(str(cooccur_row[measure]))
        elif reverse_taxa_combo in cooccur_tab.index:
            sys.exit(f"Error - Reverse taxa_combo found in co-occurrence table: {reverse_taxa_combo}.")
        else:
            cooccur_values = ["NA"] * len(measures)

        hgt_cols = '\t'.join(hgt_values)
        cooccur_cols = '\t'.join(cooccur_values)
        print(f"{taxa_combo}\t{taxon_i}\t{taxon_j}\t{tip_dist_val}\t{diff_level}\t{hgt_cols}\t{cooccur_cols}")

if __name__ == '__main__':
    main()
