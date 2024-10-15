rm(list = ls(all.names = TRUE))

source('~/scripts/ocean_mag_hgt/scripts/networks/network_functions.R')

metaG_rpkm <- read.table("/mfs/gdouglas/projects/ocean_mags/networks/combined_tables/metaG_rpkm_allsamples.tsv.gz",
                         header = TRUE, sep = "\t", row.names = 1)

run_netcomi_propr(in_tab=metaG_rpkm, random_seed=123456, outfile_gizpped="/mfs/gdouglas/projects/ocean_mags/networks/allsamples/metaG_propr_rpkm.allsamples.tsv.gz")
