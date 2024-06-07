#!/usr/bin/python3

import argparse
import sys
import pandas as pd
import os


def heatmap_tallies(tab,
                    var1,
                    var2,
                    var1_min,
                    var1_max,
                    var2_min,
                    var2_max,
                    outfile_prefix,
                    rm_var1_NA=True,
                    replace_var1_NA=False,
                    replace_var1_NA_val=0,
                    rm_var2_NA=True,
                    replace_var2_NA=False,
                    replace_var2_NA_val=0):

    outfile = outfile_prefix + "." + var1 + "." + var2 + ".tsv"

    # Check whether both settings are True or False.
    if rm_var1_NA == replace_var1_NA:
        sys.exit("rm_var1_NA and replace_var1_NA cannot both be True or False.")
    if rm_var2_NA == replace_var2_NA:
        sys.exit("rm_var2_NA and replace_var2_NA cannot both be True or False.")

    if rm_var1_NA:
        tab = tab[tab[var1].notna()]
    elif replace_var1_NA:
        tab.loc[tab[var1].isna(), var1] = replace_var1_NA_val

    if rm_var2_NA:
        tab = tab[tab[var2].notna()]
    elif replace_var2_NA:
        tab.loc[tab[var2].isna(), var1] = replace_var2_NA_val

    # Confirm that both columns are numeric.
    if not pd.api.types.is_numeric_dtype(tab[var1]):
        sys.exit(var1 + " is not numeric.")
    if not pd.api.types.is_numeric_dtype(tab[var2]):
        sys.exit(var2 + " is not numeric.")

    # Check that no rows are out side the min and max values.
    if tab[(tab[var1] < var1_min) | (tab[var1] > var1_max)].shape[0] > 0:
        sys.exit("Rows outside of " + var1 + " min and max values.")
    if tab[(tab[var2] < var2_min) | (tab[var2] > var2_max)].shape[0] > 0:
        print(tab[(tab[var2] < var2_min) | (tab[var2] > var2_max)])
        sys.exit("Rows outside of " + var2 + " min and max values.")

    var1_dist_window = (var1_max - var1_min) / 100
    var1_lower_bound = var1_min

    var2_dist_window = (var2_max - var2_min) / 100

    with open(outfile, 'w') as out_fh:
        out_fh.write(var1 + "\t" + var2 + "\tCount\n")

        for i in range(101):
            var1_upper_bound = var1_lower_bound + var1_dist_window
            tab_subset = tab[(tab[var1] >= var1_lower_bound) & (tab[var1] < var1_upper_bound)]
            var1_lower_bound_str = str(var1_lower_bound)

            var2_lower_bound = var2_min
            for j in range(101):
                var2_upper_bound = var2_lower_bound + var2_dist_window

                tab_subset_j = tab_subset[(tab_subset[var2] >= var2_lower_bound) & (tab_subset[var2] < var2_upper_bound)]

                out_fh.write(var1_lower_bound_str + "\t" + str(var2_lower_bound) + "\t" + str(tab_subset_j.shape[0]) + "\n")

                var2_lower_bound = var2_upper_bound

            var1_lower_bound = var1_upper_bound


def main():

    parser = argparse.ArgumentParser(

        description='''
Parse co-occurrence table, tips distances, MAG taxa levels, and HGT summary.
Output a table (to STDOUT) with the key information from these tables for each tip pair present.
''',

        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('--combined_in', metavar="TABLE", type=str,
                        help="Path to gzipped, tab-delimited combined table of co-occurrence, HGT, and tip distances.",
                        required=True)

    parser.add_argument('--cooccur_measure', metavar="MEASURE", type=str,
                        help="Column name from co-occurrence table to use as association measure.",
                        required=True)

    # Min and max values for co-occurrence measure.
    parser.add_argument('--min_value', metavar="MIN", type=float,
                        help="Minimum value for co-occurrence measure.",
                        required=True)

    parser.add_argument('--max_value', metavar="MAX", type=float,
                        help="Maximum value for co-occurrence measure.",
                        required=True)

    parser.add_argument('--outfolder', metavar="OUTPUT", type=str,
                        help="Output folder.",
                        required=True)

    parser.add_argument('--corr_replace_NA', action='store_true',
                        help="Replace NA values of co-occurrence measure with 0.")

    args = parser.parse_args()

    tip_dist_min = 0
    tip_dist_max = 4.04

    # For simplicity, decided to set upper limit to 100 for all measures.
    # hgt_overall_max_values = {}
    # hgt_overall_max_values["95_hit_count"] = 653
    # hgt_overall_max_values["99_hit_count"] = 781
    # hgt_overall_max_values["both_hit_count"] = 849
    # hgt_overall_max_values["95_gene_count"] = 4345
    # hgt_overall_max_values["99_gene_count"] = 6709
    # hgt_overall_max_values["both_gene_count"] = 6709
    # hgt_genus_and_above_max_values = {}
    # hgt_genus_and_above_max_values["95_hit_count"] = 222
    # hgt_genus_and_above_max_values["99_hit_count"] = 471
    # hgt_genus_and_above_max_values["both_hit_count"] = 538
    # hgt_genus_and_above_max_values["95_gene_count"] = 942
    # hgt_genus_and_above_max_values["99_gene_count"] = 2730.5
    # hgt_genus_and_above_max_values["both_gene_count"] = 2831.5

    combined = pd.read_table(filepath_or_buffer=args.combined_in, sep='\t',
                             compression='gzip', header=0, index_col=0, na_values=['NA'])

    cooccur_col = args.cooccur_measure
    if cooccur_col not in combined.columns:
        sys.exit("Column " + cooccur_col + " not found in table.")

    hgt_categories = ["95_hit_count", "99_hit_count", "both_hit_count", "95_gene_count", "99_gene_count", "both_gene_count"]
    for hgt_category in hgt_categories:
        if hgt_category not in combined.columns:
            sys.exit("Column " + hgt_category + " not found in table.")
        # Set all values above 100 to 100.
        combined.loc[combined[hgt_category] > 100, hgt_category] = 100

    if not os.path.exists(args.outfolder):
        os.makedirs(args.outfolder)

    # Subset combined table to only include rows with diff_tax_level besides Strain and Species.
    combined_higher = combined[(combined["diff_tax_level"] != "Strain") & (combined["diff_tax_level"] != "Species")]

    corr_rm_NA_flag = not args.corr_replace_NA
    corr_replace_NA_flag = args.corr_replace_NA

    out_prefix_all = args.outfolder + "/heatmap_tallies_all_tax"
    out_prefix_higher = args.outfolder + "/heatmap_tallies_genus_and_above"

    heatmap_tallies(tab=combined,
                    var1="tip_dist",
                    var2=cooccur_col,
                    rm_var2_NA=corr_rm_NA_flag,
                    replace_var2_NA=corr_replace_NA_flag,
                    replace_var2_NA_val=0,
                    var1_max=tip_dist_max,
                    var1_min=tip_dist_min,
                    var2_max=args.max_value,
                    var2_min=args.min_value,
                    outfile_prefix=out_prefix_all)

    heatmap_tallies(tab=combined_higher,
                    var1="tip_dist",
                    var2=cooccur_col,
                    rm_var2_NA=corr_rm_NA_flag,
                    replace_var2_NA=corr_replace_NA_flag,
                    replace_var2_NA_val=0,
                    var1_max=tip_dist_max,
                    var1_min=tip_dist_min,
                    var2_max=args.max_value,
                    var2_min=args.min_value,
                    outfile_prefix=out_prefix_higher)

    for hgt_category in hgt_categories:

        heatmap_tallies(tab=combined,
                        var1="tip_dist",
                        var2=hgt_category,
                        var1_max=tip_dist_max,
                        var1_min=tip_dist_min,
                        var2_max=100,
                        var2_min=0,
                        outfile_prefix=out_prefix_all)

        heatmap_tallies(tab=combined,
                        var1=cooccur_col,
                        var2=hgt_category,
                        rm_var1_NA=corr_rm_NA_flag,
                        replace_var1_NA=corr_replace_NA_flag,
                        replace_var1_NA_val=0,
                        var1_max=args.max_value,
                        var1_min=args.min_value,
                        var2_max=100,
                        var2_min=0,
                        outfile_prefix=out_prefix_all)

        heatmap_tallies(tab=combined_higher,
                        var1="tip_dist",
                        var2=hgt_category,
                        var1_max=tip_dist_max,
                        var1_min=tip_dist_min,
                        var2_max=100,
                        var2_min=0,
                        outfile_prefix=out_prefix_higher)

        heatmap_tallies(tab=combined_higher,
                        var1=cooccur_col,
                        var2=hgt_category,
                        rm_var1_NA=corr_rm_NA_flag,
                        replace_var1_NA=corr_replace_NA_flag,
                        replace_var1_NA_val=0,
                        var1_max=args.max_value,
                        var1_min=args.min_value,
                        var2_max=100,
                        var2_min=0,
                        outfile_prefix=out_prefix_higher)


if __name__ == '__main__':
    main()
