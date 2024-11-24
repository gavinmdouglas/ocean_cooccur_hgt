rm(list = ls(all.names = TRUE))

# Quick plot showing the (overall) summary of COG category expectations
# based on results that were very helpfully collated in Dmitrijeva et al. 2024.

library(ggplot2)

expected_hgt <- read.table('/mfs/gdouglas/projects/ocean_hgt_zenodo/mapfiles/Dmitrijeva2024_COG_category_HGT_summary.tsv.gz',
                           header = TRUE, sep = '\t', stringsAsFactors = FALSE)

expected_hgt <- expected_hgt[order(expected_hgt$Custom_summary), ]
