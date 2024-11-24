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

    parser.add_argument('--cooccur_measure', metavar="MEASURE", type=str,
                        help="Column name from co-occurrence table to use as association measure.",
                        required=True)

    parser.add_argument('--outprefix', metavar="PREFIX", type=str,
                        help="Output prefix.",
                        required=True)

    parser.add_argument('--hgt_measure', metavar="MEASURE", type=str,
                        help="Column name from co-occurrence table to use as association measure.",
                        required=False,
                        default="both_gene_count")

    parser.add_argument('--corr_replace_NA', action='store_true',
                        help="Replace NA values of co-occurrence measure with 0.")

    parser.add_argument('--include_lower_taxa', action='store_true',
                        help="Option to include Species and Strain levels for correlations (otherwise will be removed).")

    parser.add_argument('--extra_col', metavar="COLNAME", type=str,
                        help="Extra column name to include in output. Must have only a single value per taxon.",
                        required=False,
                        default=None)

    parser.add_argument('-m', '--min_taxa', metavar="MIN_TAXA", type=int,
                        help="Minimum number of HGT partners to consider taxon.",
                        required=False,
                        default=10)

    args = parser.parse_args()

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

    pairwise_spearman_raw = []
    partial_spearman_raw = []

    for taxon in unique_taxa:
        taxon_rows = combined[(combined["taxon_i"] == taxon) | (combined["taxon_j"] == taxon)]
        if taxon_rows.shape[0] < args.min_taxa:
            continue

        # Check that each column has at least two unique values.
        if taxon_rows[cooccur].nunique() < 2 or taxon_rows[hgt].nunique() < 2 or taxon_rows["tip_dist"].nunique() < 2:
            print("Skipping taxon " + taxon + " due to insufficient unique values.", file=sys.stderr)
            continue

        spearman_corr = pg.pairwise_corr(taxon_rows, columns=[cooccur, hgt, "tip_dist"], method='spearman')
        spearman_corr_partial = pg.partial_corr(taxon_rows, x=cooccur, y=hgt, covar="tip_dist", method='spearman')

        # Add taxon to each dataframe
        spearman_corr["taxon"] = taxon
        spearman_corr_partial["taxon"] = taxon

        if args.extra_col is not None:
            extra_col_val = taxon_rows[args.extra_col].unique()
            if len(extra_col_val) != 1:
                print("Error: Multiple values found for " + args.extra_col + " for taxon " + taxon, file=sys.stderr)
                sys.exit(1)
            extra_col_val = extra_col_val[0]
            spearman_corr[args.extra_col] = extra_col_val
            spearman_corr_partial[args.extra_col] = extra_col_val

        pairwise_spearman_raw.append(spearman_corr)
        partial_spearman_raw.append(spearman_corr_partial)

    # Combine all dataframes into one.
    pairwise_spearman = pd.concat(pairwise_spearman_raw)
    partial_spearman = pd.concat(partial_spearman_raw)

    pairwise_spearman.to_csv(args.outprefix + ".pairwise.tsv", sep='\t', index=False)
    partial_spearman.to_csv(args.outprefix + ".partial.tsv", sep='\t', index=False)


if __name__ == '__main__':
    main()
