rm(list = ls(all.names = TRUE))

source('~/scripts/ocean_mag_hgt/scripts/networks/network_functions.R')

# Run for Tara and Geotraces datasets.

tara_presence <- read.table(file = "/mfs/gdouglas/projects/ocean_mags/networks/combined_tables/metaG_presence_tara.tsv.gz",
                             header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)

tara_presence_NULL <- hyperg_cooccur_parallel(in_df = tara_presence,
                                             output_gzip_file = "/mfs/gdouglas/projects/ocean_mags/networks/dataset_subset/metaG_hyperg_cooccur.Tara.tsv.gz",
                                             num_cores = 64)


geotraces_presence <- read.table(file = "/mfs/gdouglas/projects/ocean_mags/networks/combined_tables/metaG_presence_geotraces.tsv.gz",
                            header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)

geotraces_presence_NULL <- hyperg_cooccur_parallel(in_df = geotraces_presence,
                                              output_gzip_file = "/mfs/gdouglas/projects/ocean_mags/networks/dataset_subset/metaG_hyperg_cooccur.geotraces.tsv.gz",
                                              num_cores = 64)
