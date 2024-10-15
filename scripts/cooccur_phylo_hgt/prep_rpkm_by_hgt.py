#!/usr/bin/python3

import argparse
import pandas as pd
from collections import defaultdict, Counter
import gzip
import sys


def main():

    parser = argparse.ArgumentParser(

        description='''
Prep contingency tables for COG categories enriched among BLAST hits between co-occurring and non-cooccuring genomes.
By taxonomic level and identity cut-off separately.
''',

        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--hgt_tab", metavar="HGT_TAB", type=str,
                        help="Path to HGT pairwise tallies (gzipped).",
                        required=False,
                        default="/mfs/gdouglas/projects/ocean_mags/blast_output/pairwise_tally_summary.tsv.gz")

    parser.add_argument('--rpkm_tab', metavar="MEASURE", type=str,
                        help="Gzipped coverm table of RPKM values.",
                        required=False,
                        default="/mfs/gdouglas/projects/ocean_mags/coverm/combined_tables/metaG_rpkm_allsamples.tsv.gz")

    parser.add_argument('--presence_tab', metavar="MEASURE", type=str,
                        help="Gzipped coverm table of presence/absence.",
                        required=False,
                        default="/mfs/gdouglas/projects/ocean_mags/coverm/combined_tables/metaG_presence_allsamples.tsv.gz")

    args = parser.parse_args()

    # Get all column means of table (but ignoring all 0 values).
    rpkm = pd.read_table(filepath_or_buffer=args.rpkm_tab, compression='gzip', header=0, index_col=0)
    rpkm_means = rpkm.apply(lambda x: x[x > 0].mean())

    presence = pd.read_table(filepath_or_buffer=args.presence_tab, compression='gzip', header=0, index_col=0)
    prevalence = presence.apply(lambda x: x[x > 0].count())

    rpkm_present = rpkm.where(presence > 0, 0)
    rpkm_present_means = rpkm_present.apply(lambda x: x[x > 0].mean())

    hgt_tab = pd.read_table(filepath_or_buffer=args.hgt_tab, compression='gzip', header=0, index_col=None, na_values=['NA'])
    hgt_tab = hgt_tab[hgt_tab["Highest_tax_diff"] != "Strain"]
    hgt_tab = hgt_tab[hgt_tab["Highest_tax_diff"] != "Species"]
    hgt_tab = hgt_tab.drop(columns=["95_hit_count", "99_hit_count", "both_hit_count"])

    # Split Taxa_combo column into separate columns.
    taxa1 = []
    taxa2 = []
    for index, row in hgt_tab.iterrows():
        taxa_split = row["Taxa_combo"].split(',')
        taxa1.append(taxa_split[0])
        taxa2.append(taxa_split[1])
    hgt_tab["taxon1"] = taxa1
    hgt_tab["taxon2"] = taxa2

    # Make sure both taxa are found as name in rpkm_means.
    hgt_tab = hgt_tab[hgt_tab["taxon1"].isin(rpkm_means.index)]
    hgt_tab = hgt_tab[hgt_tab["taxon2"].isin(rpkm_means.index)]

    hgt_category = ["95_gene_count", "99_gene_count", "both_gene_count"]

    higher_levels = ["Domain", "Phylum", "Class", "Order", "Family", "Genus"]

    raw_out = []

    for level in higher_levels:
        hgt_tab_level = hgt_tab[hgt_tab["Highest_tax_diff"] == level]
        for category in hgt_category:
            hgt_tab_level_category = hgt_tab_level[hgt_tab_level[category] > 0]

            all_taxa = list(hgt_tab_level_category["taxon1"]) + list(hgt_tab_level_category["taxon2"])
            taxa_counts = Counter(all_taxa)
            taxa_counts_df = pd.DataFrame(taxa_counts.items(), columns=['Taxon', 'Num_HGT_pairs'])

            taxa_counts_df["RPKM"] = taxa_counts_df["Taxon"].map(rpkm_means)
            taxa_counts_df["RPKM_present"] = taxa_counts_df["Taxon"].map(rpkm_present_means)
            taxa_counts_df["Samples_present"] = taxa_counts_df["Taxon"].map(prevalence)
            taxa_counts_df["Category"] = category
            taxa_counts_df["Level"] = level

            raw_out.append(taxa_counts_df)

    combined_out = pd.concat(raw_out, axis=0)

    combined_out.to_csv(sys.stdout, sep='\t', index=False)


if __name__ == '__main__':
    main()
