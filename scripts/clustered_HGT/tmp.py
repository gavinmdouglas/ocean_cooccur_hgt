#!/usr/bin/python3

import argparse
import os
import sys
import gzip
from collections import defaultdict
import hammingdist
import contextlib

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from functions import read_fasta

def main():

    parser = argparse.ArgumentParser(

        description='''
Parse aligned FASTAs of sequences from the same cluster.
Identify reciprocal best-hits above 95% (that are at least between genomes of different genera or above).
Also, consider all hits initially, but only output reciprocal best hits that are on passing scaffolds.
''',

formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-i', '--input',
                        metavar="FILES", type=str,
                        help="Path to file with paths to input FASTAs. One per line. Done this way to make it trivial to parallelize.",
                        required=True)

    parser.add_argument('-t', '--tax',
                        metavar="TAXTAB", type=str,
                        help="Path to (gzipped) table with taxonomic information per genome.",
                        required=False,
                        default='/mfs/gdouglas/projects/ocean_mags/mapfiles/MAG_taxa_breakdown.tsv.gz')

    parser.add_argument('-p', '--passing',
                        metavar="PASS", type=str,
                        help="Path to file with passing scaffolds.",
                        required=False,
                        default='/mfs/gdouglas/projects/ocean_mags/Sunagawa_dataset/scaffolds_5000bp.txt')

    parser.add_argument('-g', '--geneinfo',
                        metavar="GENEINFO", type=str,
                        help="Path to (gzipped) table with gene information per genome.",
                        required=False,
                        default='/mfs/gdouglas/projects/ocean_mags/Sunagawa_dataset/gene_info_allscaffolds.tsv.gz')

    args = parser.parse_args()

    genes = set()
    with open(args.input, 'r') as input_fh:
        for filepath in input_fh:
            filepath = filepath.strip()
            if not filepath or len(filepath) == 0:
                continue
            seqs = read_fasta(filepath)
            for seq in seqs:
                genes.add(seq)

    with gzip.open(args.geneinfo, "rb") as gene_fh:
        gene_header = gene_fh.readline().decode("utf-8").strip().split("\t")
        print('\t'.join(gene_header))

        for gene_line in gene_fh:
            gene_line = gene_line.decode("utf-8").strip().split("\t")
            if gene_line[0] in genes:
                print('\t'.join(gene_line))

if __name__ == '__main__':
    main()
