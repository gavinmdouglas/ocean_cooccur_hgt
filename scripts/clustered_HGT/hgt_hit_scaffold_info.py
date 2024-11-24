#!/usr/bin/python3

import argparse
import os
import sys
import gzip

gene_to_scaffold = {}

with gzip.open('/mfs/gdouglas/projects/ocean_mags/Sunagawa_dataset/gene_info_allscaffolds.tsv.gz', 'rt') as gene_info_fh:
    gene_info_fh.readline()
    for line in gene_info_fh:
        line = line.strip().split('\t')
        gene_to_scaffold[line[0]] = line[2]


with gzip.open('/mfs/gdouglas/projects/ocean_hgt_zenodo/putative_hgt/cluster/all_best_hits.tsv.gz', 'rt') as cluster_fh:
    cluster_fh.readline()
    for line in cluster_fh:
        line = line.strip().split('\t')
        gene1 = line[0]
        gene2 = line[1]
        scaffold1 = gene_to_scaffold[gene1]
        scaffold2 = gene_to_scaffold[gene2]



/mfs/gdouglas/projects/ocean_mags/Sunagawa_dataset/gene_info_allscaffolds.tsv.gz