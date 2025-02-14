#!/usr/bin/python3

import gzip
import argparse
import sys
import pandas as pd
from itertools import combinations



with gzip.open("/mfs/gdouglas/projects/ocean_mags/progenomes_analyses/networks/combined_tables/metaG_presence_allsamples.tsv.gz", 'rt') as coverm_in:
    header = coverm_in.readline().strip().split('\t')
    if header[0] != 'sample':
        sys.exit("First column of coverm table must be 'sample'")
    taxa = sorted(header[1:])
taxa_combinations = list(combinations(taxa, 2))

print(taxa_combinations)