rm(list = ls(all.names = TRUE))

source('~/scripts/ocean_mag_hgt/scripts/networks/network_functions.R')

# Run for MGS samples filtered to be enriched for free-living cells (or less filtered, so they have more particle-associated cells).
metaG_rpkm_freeliv <- read.table(file = "/mfs/gdouglas/projects/ocean_mags/networks/combined_tables/metaG_rpkm_freeliv.tsv.gz",
                                     header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)

metaG_rpkm_lessfiltered <- read.table(file = "/mfs/gdouglas/projects/ocean_mags/networks/combined_tables/metaG_rpkm_lessfiltered.tsv.gz",
                                          header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)

run_netcomi_propr(in_tab=metaG_rpkm_freeliv, random_seed=41, outfile_gizpped="/mfs/gdouglas/projects/ocean_mags/networks/filtersplit/metaG_propr_rpkm.freeliv.tsv.gz")

run_netcomi_propr(in_tab=metaG_rpkm_lessfiltered, random_seed=51, outfile_gizpped="/mfs/gdouglas/projects/ocean_mags/networks/filtersplit/metaG_propr_rpkm.lessfiltered.tsv.gz")
