#!/usr/bin/python3

import argparse
import sys
import pandas as pd
import os


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

    parser.add_argument('--hgt_measure', metavar="MEASURE", type=str,
                        help="Column name from co-occurrence table to use as association measure.",
                        required=False,
                        default="both_gene_count")

    args = parser.parse_args()

    combined = pd.read_table(filepath_or_buffer=args.combined_in, sep='\t',
                             compression='gzip', header=0, index_col=0, na_values=['NA'])

    cooccur = args.cooccurrence_measure
    hgt = args.hgt_measure

    combined = combined[["taxon_i", "taxon_j", "tip_dist", cooccur, hgt]]

    combined = combined[combined["tip_dist"].notna()]
    combined = combined[combined[cooccur].notna()]
    combined = combined[combined[hgt] >= 1]

    


if __name__ == '__main__':
    main()
