#!/usr/bin/python3

import argparse
import sys
import pandas as pd

filepath = "/mfs/gdouglas/projects/ocean_mags/coverm/network_working_allsamples/metaG_simple_cooccur_prepped_tabs_with_env_mean_diffs/Taxa_4934.tsv"

# Increase pandas print width.
pd.set_option('display.width', 1000)
pd.set_option('display.max_columns', None)

test = pd.read_table(filepath_or_buffer=filepath, sep='\t', header=0, index_col=0)
print(test)
test["cooccur_simple_cooccur"] = test["cooccur_simple_cooccur"].sample(frac=1).values

print(test)


