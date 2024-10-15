#!/usr/bin/python3

import argparse
import sys
import pandas as pd
import pingouin as pg


def main():

    parser = argparse.ArgumentParser(

        description='''
Perform Spearman correlation between all pairwise variables of interest.
''',

        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('--combined_in', metavar="TABLE", type=str,
                        help="Path to gzipped, tab-delimited combined table of co-occurrence, HGT, and tip distances.",
                        required=True)

    parser.add_argument('--mean_env_diff', metavar="TABLE", type=str,
                        help="Path to gzipped, tab-delimited table of pairwise differences between mean env. measures between genomes.",
                        required=True)

    parser.add_argument('--cooccur_measure', metavar="MEASURE", type=str,
                        help="Column name from co-occurrence table to use as association measure.",
                        required=True)

    parser.add_argument('--output', metavar="FOLDER", type=str,
                        help="Output folder.",
                        required=True)

    parser.add_argument('--hgt_measure', metavar="MEASURE", type=str,
                        help="Column name from co-occurrence table to use as association measure.",
                        required=False,
                        default="both_gene_count")

    parser.add_argument('--include_lower_taxa', action='store_true',
                        help="Option to include Species and Strain levels for correlations (otherwise will be removed).")

    parser.add_argument('--corr_replace_NA', action='store_true',
                        help="Replace NA values of co-occurrence measure with 0.")

    parser.add_argument('--extra_col', metavar="COLNAME", type=str,
                        help="Extra column name to include in output. Must have only a single value per taxon.",
                        required=False,
                        default=None)

    args = parser.parse_args()

    mean_env_diff = pd.read_table(filepath_or_buffer=args.mean_env_diff, sep='\t',
                                  compression='gzip', header=0, index_col=0, na_values=['NA'])
    # Check if any rows have NA values.
    if mean_env_diff.isna().sum().sum() > 0:
        print("Warning: NA values found in mean environment difference table. Removing these rows.", 
                file=sys.stderr)

    mean_env_diff = mean_env_diff.dropna()

    combined = pd.read_table(filepath_or_buffer=args.combined_in, sep='\t',
                             compression='gzip', header=0, index_col=0, na_values=['NA'])

    cooccur = args.cooccur_measure
    hgt = args.hgt_measure

    if not args.include_lower_taxa:
        # Remove where diff_tax_level is "Strain" or "Species".
        combined = combined[combined["diff_tax_level"] != "Strain"]
        combined = combined[combined["diff_tax_level"] != "Species"]

    # Perform additional filtering.
    if args.extra_col is not None:
        if args.extra_col not in combined.columns:
            print("Error: " + args.extra_col + " not found in combined table.", file=sys.stderr)
            sys.exit(1)
        combined = combined[["taxon_i", "taxon_j", "tip_dist", cooccur, hgt, args.extra_col]]
    else:
        combined = combined[["taxon_i", "taxon_j", "tip_dist", cooccur, hgt]]

    combined = combined[combined["tip_dist"].notna()]
    combined = combined[combined[hgt].notna()]
    combined = combined[combined[hgt] >= 1]

    num_cooccur_na = combined[cooccur].isna().sum()
    if num_cooccur_na > 0:
        if not args.corr_replace_NA:
            print(str(num_cooccur_na) + " NA values found in co-occurrence measure. Removing these rows.", file=sys.stderr)
            combined = combined[combined[cooccur].notna()]
        else:
            print(str(num_cooccur_na) + " NA values found in co-occurrence measure. Filling in with 0.")
            combined.fillna({cooccur: 0}, inplace=True)


    # Get all unique values of combined taxon_i and taxon_j columns.
    taxon_i_values = combined["taxon_i"].unique()
    taxon_j_values = combined["taxon_j"].unique()
    unique_taxa = list(set(taxon_i_values) | set(taxon_j_values))

    for taxon in unique_taxa:
        taxon_rows = combined[(combined["taxon_i"] == taxon) | (combined["taxon_j"] == taxon)]

        subset_taxa_combos = taxon_rows.index

        # Figure out which of these combos is in the mean_env_diff table.
        mean_env_diff_combos = mean_env_diff.index
        intersecting_combos = [combo for combo in subset_taxa_combos if combo in mean_env_diff_combos]

        if len(taxon_rows[hgt].unique()) < 5:
            continue

        if len(intersecting_combos) >= 10:
            mean_env_diff_values = mean_env_diff.loc[intersecting_combos]
            combined_values = combined.loc[intersecting_combos]

            # Merge tables.
            merged = pd.merge(left=mean_env_diff_values, right=combined_values, left_index=True, right_index=True)

            # Write out table.
            merged.to_csv(path_or_buf=args.output + "/" + taxon + ".tsv",
                          sep='\t', index=True, header=True)


if __name__ == '__main__':
    main()
