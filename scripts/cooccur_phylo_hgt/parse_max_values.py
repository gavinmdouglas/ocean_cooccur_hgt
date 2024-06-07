#!/usr/bin/python3

import gzip
import pandas as pd

tip_dist_file = "/mfs/gdouglas/projects/ocean_mags/phylogenetic_analyses/tip_dist.tsv.gz"
hgt_file = "/mfs/gdouglas/projects/ocean_mags/blast_output/pairwise_tally_summary.tsv.gz"

# Read in tip distances and get max value in table.
tip_dist = pd.read_csv(tip_dist_file, sep="\t", index_col=0)
tip_dist_max = tip_dist.max().max()

del tip_dist

print("measure\tmin\tmax")
print("tip_dist\t0\t" + str(tip_dist_max))

# Get max value in HGT file.
hgt_categories = ["95_hit_count", "99_hit_count", "both_hit_count", "95_gene_count", "99_gene_count", "both_gene_count"] 
hgt_max_overall = {}
hgt_max_genus_and_above = {}

for hgt_category in hgt_categories:
    hgt_max_overall[hgt_category] = -1
    hgt_max_genus_and_above[hgt_category] = -1

with gzip.open(hgt_file, "rt") as hgt_fh:
    hgt_header = hgt_fh.readline().strip().split("\t")
    hgt_col_to_i = {col: i for i, col in enumerate(hgt_header)}
    for hgt_line in hgt_fh:
        hgt_line = hgt_line.strip().split("\t")

        for hgt_category in hgt_categories:
            if hgt_line[hgt_col_to_i[hgt_category]] == "NA":
                    continue
            hgt_max_overall[hgt_category] = max(hgt_max_overall[hgt_category], float(hgt_line[hgt_col_to_i[hgt_category]]))
        
        if hgt_line[hgt_col_to_i["Highest_tax_diff"]] != "Species" and hgt_line[hgt_col_to_i["Highest_tax_diff"]] != "Strain":
            for hgt_category in hgt_categories:
                if hgt_line[hgt_col_to_i[hgt_category]] == "NA":
                    continue
                hgt_max_genus_and_above[hgt_category] = max(hgt_max_genus_and_above[hgt_category], float(hgt_line[hgt_col_to_i[hgt_category]]))

for hgt_category in hgt_categories:
    print("hgt_overall|" + hgt_category + "\t0\t" + str(hgt_max_overall[hgt_category]))
    print("hgt_genus_and_above|" + hgt_category + "\t0\t" + str(hgt_max_genus_and_above[hgt_category]))

# Get min and max values in network files.
network_dir = "/mfs/gdouglas/projects/ocean_mags/coverm/network_working/"

hyperg_max = -1
with gzip.open(network_dir + "metaG_hyperg_cooccur.tsv.gz", "rt") as hyperg_fh:
    hyperg_header = hyperg_fh.readline().strip().split("\t")
    hyperg_col_to_i = {col: i for i, col in enumerate(hyperg_header)}
    for hyperg_line in hyperg_fh:
        hyperg_line = hyperg_line.strip().split("\t")
        if hyperg_line[6] == "NA":
            continue
        hyperg_max = max(hyperg_max, float(hyperg_line[6]))
print("metaG_hyperg\t0\t" + str(hyperg_max))

hyperg_max = -1
with gzip.open(network_dir + "metaT_hyperg_cooccur.tsv.gz", "rt") as hyperg_fh:
    hyperg_fh.readline()
    for hyperg_line in hyperg_fh:
        hyperg_line = hyperg_line.strip().split("\t")
        if hyperg_line[6] == "NA":
            continue
        hyperg_max = max(hyperg_max, float(hyperg_line[6]))
print("metaT_hyperg\t0\t" + str(hyperg_max))

# Simple co-occurrence trivially varies from 0 to 1.
print("metaG_simple\t0\t1.0")
print("metaT_simple\t0\t1.0")

propr_min = 1
propr_max = -1
with gzip.open(network_dir + "metaG_propr_rpkm.tsv.gz", "rt") as propr_fh:
    propr_header = propr_fh.readline()
    for propr_line in propr_fh:
        propr_line = propr_line.strip().split("\t")
        if propr_line[2] == "NA":
            continue
        propr_min = min(propr_min, float(propr_line[2]))
        propr_max = max(propr_max, float(propr_line[2]))
print("metaG_propr\t" + str(propr_min) + "\t" + str(propr_max))

propr_min = 1
propr_max = -1
with gzip.open(network_dir + "metaT_propr_rpkm.tsv.gz", "rt") as propr_fh:
    propr_fh.readline()
    for propr_line in propr_fh:
        propr_line = propr_line.strip().split("\t")
        if propr_line[2] == "NA":
            continue
        propr_min = min(propr_min, float(propr_line[2]))
        propr_max = max(propr_max, float(propr_line[2]))
print("metaT_propr\t" + str(propr_min) + "\t" + str(propr_max))

spieceasi_min = 1
spieceasi_max = -1
with gzip.open(network_dir + "metaG_spieceasi_rpkm.tsv.gz", "rt") as spieceasi_fh:
    spieceasi_header = spieceasi_fh.readline()
    for spieceasi_line in spieceasi_fh:
        spieceasi_line = spieceasi_line.strip().split("\t")
        if spieceasi_line[2] == "NA":
            continue
        spieceasi_min = min(spieceasi_min, float(spieceasi_line[2]))
        spieceasi_max = max(spieceasi_max, float(spieceasi_line[2]))
print("metaG_spieceasi\t" + str(spieceasi_min) + "\t" + str(spieceasi_max))

spieceasi_min = 1
spieceasi_max = -1
with gzip.open(network_dir + "metaT_spieceasi_rpkm.tsv.gz", "rt") as spieceasi_fh:
    spieceasi_fh.readline()
    for spieceasi_line in spieceasi_fh:
        spieceasi_line = spieceasi_line.strip().split("\t")
        if spieceasi_line[2] == "NA":
            continue
        spieceasi_min = min(spieceasi_min, float(spieceasi_line[2]))
        spieceasi_max = max(spieceasi_max, float(spieceasi_line[2]))
print("metaT_spieceasi\t" + str(spieceasi_min) + "\t" + str(spieceasi_max))
