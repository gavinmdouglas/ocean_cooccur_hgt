rm(list = ls(all.names = TRUE))

source('~/scripts/ocean_mag_hgt/scripts/networks/network_functions.R')

metaG_presence <- read.table(file = "/mfs/gdouglas/projects/ocean_mags/networks/combined_tables/metaG_presence_allsamples.tsv.gz",
                             header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)

metaG_hyperg_NULL <- hyperg_cooccur_parallel(in_df = metaG_presence,
                                             output_gzip_file = "/mfs/gdouglas/projects/ocean_mags/networks/allsamples/metaG_hyperg_cooccur.allsamples.tsv.gz",
                                             num_cores = 64)
