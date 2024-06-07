#!/usr/bin/python3

import argparse
from itertools import combinations
import pandas as pd
import numpy as np


def main():

    parser = argparse.ArgumentParser(

        description='''
Compute pairwise mean (absolute) differences between genomes based on mean environmental data.
''',

        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('--tab', metavar="IN_TABLE", type=str,
                        help="Path to (gzipped) tab-delimited table of mean environmental data per genome.",
                        required=True)

    args = parser.parse_args()

    tab = pd.read_table(filepath_or_buffer=args.tab, compression='gzip', header=0, index_col=0, na_values="NA")

    taxa_combinations = list(combinations(sorted(list(tab.index)), 2))

    print("\t".join(["taxa_combo", "taxon_i", "taxon_j"] + list(tab.columns)))

    for i in range(len(taxa_combinations)):
        taxon_i = taxa_combinations[i][0]
        taxon_j = taxa_combinations[i][1]
        outline = [",".join(sorted([taxon_i, taxon_j])), taxon_i, taxon_j]
        for col in tab.columns:
            if np.isnan(tab.loc[taxon_i, col]) or np.isnan(tab.loc[taxon_j, col]):
                outline.append("NA")
            else:
                outline.append(str(np.abs(tab.loc[taxon_i, col] - tab.loc[taxon_j, col])))

        print("\t".join(outline))

if __name__ == '__main__':
    main()
