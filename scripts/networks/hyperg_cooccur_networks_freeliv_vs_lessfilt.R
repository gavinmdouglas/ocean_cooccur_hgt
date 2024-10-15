rm(list = ls(all.names = TRUE))

source('~/scripts/ocean_mag_hgt/scripts/networks/network_functions.R')

# Run for MGS samples filtered to be enriched for free-living cells (or less filtered, so they have more particle-associated cells).

metaG_presence_freeliv <- read.table(file = "/mfs/gdouglas/projects/ocean_mags/networks/combined_tables/metaG_presence_freeliv.tsv.gz",
                                     header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)

metaG_presence_freeliv_hyperg_NULL <- hyperg_cooccur_parallel(in_df = metaG_presence_freeliv,
                                                              output_gzip_file = "/mfs/gdouglas/projects/ocean_mags/networks/filtersplit/metaG_hyperg_cooccur.freeliv.tsv.gz",
                                                              num_cores = 64)

metaG_presence_lessfiltered <- read.table(file = "/mfs/gdouglas/projects/ocean_mags/networks/combined_tables/metaG_presence_lessfiltered.tsv.gz",
                                     header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)

metaG_presence_lessfiltered_hyperg_NULL <- hyperg_cooccur_parallel(in_df = metaG_presence_lessfiltered,
                                                              output_gzip_file = "/mfs/gdouglas/projects/ocean_mags/networks/filtersplit/metaG_hyperg_cooccur.lessfiltered.tsv.gz",
                                                              num_cores = 64)
