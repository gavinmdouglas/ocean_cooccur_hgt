#!/usr/bin/python3

import gzip
import argparse
import sys
import pandas as pd

def main():

    parser = argparse.ArgumentParser(

        description='''
Parse co-occurrence table and output a matching table with the tip distances
between each pair of taxa, and the highest taxonomic level where they differ.
Output will be written to standard output.
''',

        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-t', '--table', metavar='COOCCUR_TAB', type=str,
                        help="Path to gzipped co-occurrence table. Must be tab. delimited and start with the columns taxon_i and taxon_j",
                        required=True)

    parser.add_argument("--tip_distances", metavar="TIP_DIST", type=str,
                        help="Path to tip distances table (gzipped).",
                        required=False,
                        default="/mfs/gdouglas/projects/water_mags/water_mag_analysis/phylogenetic_analyses/tip_dist.tsv.gz")

    parser.add_argument("--taxa_tab", metavar="TAX_TAB", type=str,
                        help="Path to taxonomic table breakdown (gzipped).",
                        required=False,
                        default="/mfs/gdouglas/projects/water_mags/water_mag_analysis/mapfiles/MAG_taxa_breakdown.tsv.gz")

    parser.add_argument('--assoc', action='store_true',
                        help="Print association value from RPKM association table.")

    args = parser.parse_args()

    # Read in taxa table as pandas dataframe, and set second column as index.
    taxa_tab = pd.read_table(filepath_or_buffer=args.taxa_tab, compression='gzip', header=0, index_col=1)

    # Read in tip distances table as pandas dataframe.
    tip_dist = pd.read_table(filepath_or_buffer=args.tip_distances, compression='gzip', header=0, index_col=0)
    
    if args.assoc:
        print("taxon_i\ttaxon_j\tassoc\ttip_dist\tdiff_tax_level")
    else:
        print("taxon_i\ttaxon_j\ttip_dist\tdiff_tax_level")

    with gzip.open(args.table, 'rt') as f:
        header = f.readline().strip().split('\t')
        if (header[0] != 'taxon_i' or header[1] != 'taxon_j') and (header[0] != 'taxoni' or header[1] != 'taxonj'):
            sys.exit("First two columns must be 'taxon_i' and 'taxon_j' or 'taxoni' and 'taxonj'")

        if args.assoc:
            if "asso" not in header:
                sys.exit("Association column not found in co-occurrence table.")
            assoc_i = header.index("asso")

        for line in f:
            line = line.strip().split('\t')
            taxon_i = line[0]
            taxon_j = line[1]

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
            
            if args.assoc:
                print(f"{taxon_i}\t{taxon_j}\t{line[assoc_i]}\t{tip_dist_val}\t{diff_level}")
            else:
                print(f"{taxon_i}\t{taxon_j}\t{tip_dist_val}\t{diff_level}")


if __name__ == '__main__':
    main()

