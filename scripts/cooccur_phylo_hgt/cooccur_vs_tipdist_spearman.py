#!/usr/bin/python3

import argparse
import scipy.stats
import pandas as pd

def run_spearmanr(tab, x_col, y_col):
    spearman_out = scipy.stats.spearmanr(tab[x_col], tab[y_col])
    return (spearman_out.correlation, spearman_out.pvalue)


def main():

    parser = argparse.ArgumentParser(

        description=
        '''
        Compute Spearman correlation between tip distance and each co-occurrence metric,
        ''',

        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('--input', metavar="PATH", type=str,
                        help="Input table with all variables of interest.",
                        required=False,
                        default="/mfs/gdouglas/projects/ocean_mags/networks/allsamples/clusterbased/cooccur_and_tipdist.combined.tsv.ORIG.gz")
    
    parser.add_argument('--output', metavar="PATH", type=str,
                        help="Output table with Spearman correlations.",
                        required=False,
                        default="/mfs/gdouglas/projects/ocean_mags/networks/allsamples/clusterbased/cooccur_and_tipdist.combined.spearman.tsv")

    args = parser.parse_args()

    combined_tab = pd.read_table(filepath_or_buffer=args.input, compression='gzip', header=0, index_col=0)

    # Compute Spearman correlation between tip distance and each co-occurrence metric,
    # based on all taxa, as well as based on all taxa that differ at species and lower,
    # as well as taxa that differ at genus and above.
    tab_name = []
    cooccur_name = []
    spearman_r_values = []
    spearman_p_values = []

    prepped_tables = {}
    prepped_tables['All_levels'] = combined_tab
    prepped_tables['Species_and_lower'] = combined_tab[combined_tab['diff_tax_level'].isin(['Species', 'Strain'])]
    prepped_tables['Genus_and_above'] = combined_tab[~combined_tab['diff_tax_level'].isin(['Species', 'Strain'])]

    for prepped_tabname in sorted(prepped_tables.keys()):
        for col in ['simple_cooccur', 'hyperg_ratio', 'propr_asso', 'propr_diss', 'propr_adja']:
            prepped_tab = prepped_tables[prepped_tabname].copy()
            prepped_tab = prepped_tab.dropna(subset=[col])
            tab_name.append(prepped_tabname)
            cooccur_name.append(col)
            spearman_r, spearman_p = run_spearmanr(prepped_tab, 'tip_dist', col)
            spearman_r_values.append(spearman_r)
            spearman_p_values.append(spearman_p)

    spearman_tab = pd.DataFrame({'tab_name': tab_name,
                                 'cooccur_name': cooccur_name,
                                 'spearman_r': spearman_r_values,
                                 'spearman_p': spearman_p_values})

    spearman_tab.to_csv(args.output, sep='\t', index=False)


if __name__ == '__main__':
    main()
