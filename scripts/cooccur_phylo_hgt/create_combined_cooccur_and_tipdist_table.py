#!/usr/bin/python3

import argparse
import pandas as pd
import numpy as np
import sys


def np_close_check(a, b, tolerance=1e-6):
    if not np.isclose(a, b, rtol=tolerance).all():
        for i in range(len(a)):
            if not np.isclose(a[i], b[i], rtol=tolerance):
                print(a[i], b[i])
        sys.exit("Error - np_close_check failed.")


def main():

    parser = argparse.ArgumentParser(

        description='''
        Parse co-occurrence table, tips distances, MAG taxa levels, and HGT summary.
        Put them all into a single table (with rows with NA tips removed!)
        ''',

        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('--simple', metavar="PATH", type=str,
                        help="Prepped gzipped table with simple co-occurrence values (along with other info, including phylogenetic dist).",
                        required=False,
                        default="/mfs/gdouglas/projects/ocean_mags/networks/allsamples/clusterbased/metaG_simple_cooccur.allsamples.combined.tsv.gz")

    parser.add_argument('--hyperg', metavar="MEASURE", type=str,
                        help="Prepped gzipped table with hyperg co-occurrence values (along with other info, including phylogenetic dist).",
                        required=False,
                        default="/mfs/gdouglas/projects/ocean_mags/networks/allsamples/clusterbased/metaG_hyperg_cooccur.allsamples.combined.tsv.gz")

    parser.add_argument('--propr', metavar="MEASURE", type=str,
                        help="Prepped gzipped table with hyperg co-occurrence values (along with other info, including phylogenetic dist).",
                        required=False,
                        default="/mfs/gdouglas/projects/ocean_mags/networks/allsamples/clusterbased/metaG_propr_rpkm.allsamples.combined.tsv.gz")

    parser.add_argument('--output', metavar="PATH", type=str,
                        help="Output gzipped table.",
                        required=False,
                        default="/mfs/gdouglas/projects/ocean_mags/networks/allsamples/clusterbased/cooccur_and_tipdist.combined.tsv.gz")

    parser.add_argument('--taxa_combos_no_tipdist', metavar="PATH", type=str,
                        help="Output table with taxa combinations with NA tip distances.",
                        required=False,
                        default="/mfs/gdouglas/projects/ocean_hgt_zenodo/taxa_subsets/taxa_pairs_NA_tipdist.tsv.gz")


    args = parser.parse_args()

    simple_tab = pd.read_table(filepath_or_buffer=args.simple, compression='gzip', header=0, index_col=0)

    # Get rows with NA values in "tip_dist" column, and write out for later sanity checks of these taxa.
    simple_tab_na_tipdist = simple_tab[simple_tab['tip_dist'].isna()]
    simple_tab_na_tipdist[['taxon_i', 'taxon_j']].to_csv(args.taxa_combos_no_tipdist, sep='\t', compression='gzip', na_rep='NA')

    simple_tab = simple_tab.dropna(subset=['tip_dist'])

    # Keep only columns "tip_dist", "diff_tax_level", and "cooccur_simple_cooccur".
    simple_tab = simple_tab[['tip_dist', 'diff_tax_level', 'cooccur_simple_cooccur']]
    simple_tab = simple_tab.rename(columns={'cooccur_simple_cooccur': 'simple_cooccur'})

    hyperg_tab = pd.read_table(filepath_or_buffer=args.hyperg, compression='gzip', header=0, index_col=0)
    hyperg_tab = hyperg_tab.dropna(subset=['tip_dist'])
    hyperg_tab['cooccur_ratio'] = hyperg_tab['cooccur_obs'] / hyperg_tab['cooccur_exp']
    hyperg_tab = hyperg_tab.rename(columns={'cooccur_ratio': 'hyperg_ratio'})

    propr_tab = pd.read_table(filepath_or_buffer=args.propr, compression='gzip', header=0, index_col=0)
    propr_tab = propr_tab.dropna(subset=['tip_dist'])

    # Add "propr" to column names: cooccur_asso    cooccur_diss    cooccur_adja
    propr_tab = propr_tab.rename(columns={'cooccur_asso': 'propr_asso',
                                          'cooccur_diss': 'propr_diss',
                                          'cooccur_adja': 'propr_adja'})

    # Identify row names that match across all three tables.
    common_indices = simple_tab.index.intersection(hyperg_tab.index).intersection(propr_tab.index)

    # Sanity check that tip distances are identical across all three tables for these common indices.
    # Same for "diff_tax_level"
    np_close_check(simple_tab.loc[common_indices, 'tip_dist'], hyperg_tab.loc[common_indices, 'tip_dist'])
    np_close_check(simple_tab.loc[common_indices, 'tip_dist'], propr_tab.loc[common_indices, 'tip_dist'])
    assert all(simple_tab.loc[common_indices, 'diff_tax_level'] == hyperg_tab.loc[common_indices, 'diff_tax_level'])
    assert all(simple_tab.loc[common_indices, 'diff_tax_level'] == propr_tab.loc[common_indices, 'diff_tax_level'])

    # Merge the three tables.
    combined_tab = simple_tab.loc[common_indices, :]
    combined_tab['hyperg_ratio'] = hyperg_tab.loc[common_indices, 'hyperg_ratio']
    combined_tab['propr_asso'] = propr_tab.loc[common_indices, 'propr_asso']
    combined_tab['propr_diss'] = propr_tab.loc[common_indices, 'propr_diss']
    combined_tab['propr_adja'] = propr_tab.loc[common_indices, 'propr_adja']

    # Output the combined table.
    combined_tab.to_csv(args.output, sep='\t', compression='gzip', na_rep='NA')


if __name__ == '__main__':
    main()
