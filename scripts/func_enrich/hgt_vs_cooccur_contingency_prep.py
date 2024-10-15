#!/usr/bin/python3

import argparse
import pandas as pd
from collections import defaultdict
import sys


def main():

    parser = argparse.ArgumentParser(

        description='''
Parse combined co-occurrence and HGT table.
For each taxon, output the number of times each of the following occurs:
1. Positive co-occurrence is significant and HGT is positive.
2. Positive co-occurrence is significant and HGT is negative.
3. Co-occurrence is not significant (or is negative co-occurrence) and HGT is positive.
4. Co-occurrence is not significant (or is negative co-occurrence) and HGT is negative.
''',

        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('--combined_in', metavar='COMBINED_TAB', type=str,
                        help="Path to gzipped combined table.",
                        required=True)

    parser.add_argument('--cooccur_measure', metavar="MEASURE", type=str,
                        help="Column name to use for co-occurrence measure.",
                        required=True)

    parser.add_argument("--taxa_tab", metavar="TAX_TAB", type=str,
                        help="Path to taxonomic table breakdown (gzipped).",
                        required=False,
                        default="/mfs/gdouglas/projects/ocean_mags/mapfiles/MAG_taxa_breakdown.tsv.gz")

    parser.add_argument('--cooccur_measure_p', metavar="MEASURE", type=str,
                        help="Column name to use for co-occurrence measure P-value.",
                        required=True)

    parser.add_argument('--hgt_measure', metavar="MEASURE", type=str,
                        help="Column name to use for HGT measure",
                        required=True)

    parser.add_argument('--cooccur_cutoff', metavar="MEASURE", type=float,
                        help="Estimate cut-off to use to identify a significant hit as higher or lower than expected by chance.",
                        required=False,
                        default=1.0)

    parser.add_argument('--include_lower_taxa', action='store_true',
                        help="Option to keep Species and Strain levels (otherwise will be removed).")

    args = parser.parse_args()

    # Read in mapping of genome ID to taxonomic category for each taxonomic levels.
    taxa_map = defaultdict(dict)
    taxa_tab = pd.read_table(filepath_or_buffer=args.taxa_tab, sep='\t',
                                compression='gzip', header=0, index_col=1, na_values=['NA'])
    for taxon in taxa_tab.index.values:
        for level in taxa_tab.columns.values:
            taxa_map[taxon][level] = taxa_tab.loc[taxon, level]

    combined = pd.read_table(filepath_or_buffer=args.combined_in, sep='\t',
                             compression='gzip', header=0, index_col=0, na_values=['NA'])

    cooccur = args.cooccur_measure
    hgt = args.hgt_measure
    cooccur_p = args.cooccur_measure_p
    cutoff = args.cooccur_cutoff

    if not args.include_lower_taxa:
        # Remove where diff_tax_level is "Strain" or "Species".
        combined = combined[combined["diff_tax_level"] != "Strain"]
        combined = combined[combined["diff_tax_level"] != "Species"]

    tax_levels = combined["diff_tax_level"].unique()

    # all_unique_taxa = set(combined["taxon_i"].unique().tolist() + combined["taxon_j"].unique().tolist())
    # all_unique_taxa = sorted(list(all_unique_taxa))
    # all_taxa_combinations = []
    # all_taxa_combinations_raw = list(combinations(all_unique_taxa, 2))
    # for i, (taxon_i, taxon_j) in enumerate(all_taxa_combinations_raw):
    #     all_taxa_combinations.append(f"{taxon_i},{taxon_j}")
    # assert all_taxa_combinations == sorted(list(combined.index)), \
    # "all_taxa_combinations does not match the sorted combined index"

    print("Tax_Level\tFull_Taxonomy\tCooccur_HGT\tCooccur_No_HGT\tNo_Cooccur_HGT\tNo_Cooccur_No_HGT")

    for tax_level in tax_levels:

        combined_subset = combined[combined["diff_tax_level"] == tax_level]

        tallies_cooccur_hgt = defaultdict(int)
        tallies_cooccur_no_hgt = defaultdict(int)
        tallies_no_cooccur_hgt = defaultdict(int)
        tallies_no_cooccur_no_hgt = defaultdict(int)

        unique_taxa = set()

        # Loop over all rows in the subset.
        for i, row in combined_subset.iterrows():
            for taxon in [row["taxon_i"], row["taxon_j"]]:
                full_taxonomy = taxa_map[taxon][tax_level]
                unique_taxa.add(full_taxonomy)

                hgt_flag = False
                if row[hgt] > 0:
                    hgt_flag = True

                cooccur_flag = False
                if row[cooccur] > cutoff and row[cooccur_p] < 0.05:
                    cooccur_flag = True

                if not hgt_flag and not cooccur_flag:
                    tallies_no_cooccur_no_hgt[full_taxonomy] += 1
                elif not hgt_flag and cooccur_flag:
                    tallies_cooccur_no_hgt[full_taxonomy] += 1
                elif hgt_flag and not cooccur_flag:
                    tallies_no_cooccur_hgt[full_taxonomy] += 1
                elif hgt_flag and cooccur_flag:
                    tallies_cooccur_hgt[full_taxonomy] += 1
 
        for full_taxonomy in sorted(list(unique_taxa)):
            print(f"{tax_level}\t{full_taxonomy}\t{tallies_cooccur_hgt[full_taxonomy]}\t{tallies_cooccur_no_hgt[full_taxonomy]}\t{tallies_no_cooccur_hgt[full_taxonomy]}\t{tallies_no_cooccur_no_hgt[full_taxonomy]}")


if __name__ == '__main__':
    main()
