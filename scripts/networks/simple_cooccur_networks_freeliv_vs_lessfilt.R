rm(list = ls(all.names = TRUE))

source("~/scripts/ocean_mag_hgt/scripts/networks/network_functions.R")

# Run for MGS samples filtered to be enriched for free-living cells (or less filtered, so they have more particle-associated cells).
metaG_presence_freeliv <- read.table(file = "/mfs/gdouglas/projects/ocean_mags/networks/combined_tables/metaG_presence_freeliv.tsv.gz",
                                     header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)

metaG_presence_lessfiltered <- read.table(file = "/mfs/gdouglas/projects/ocean_mags/networks/combined_tables/metaG_presence_lessfiltered.tsv.gz",
                                          header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)

metaG_presence_freeliv_simple_NULL <- simple_cooccur_parallel(in_df = metaG_presence_freeliv,
                                                              ncores = 64,
                                                              tmp_dir = "/mfs/gdouglas/tmp/metaG_presence_freeliv_simple_cooccur_tmp/",
                                                              output_gzip_file = "/mfs/gdouglas/projects/ocean_mags/networks/filtersplit/metaG_simple_cooccur.freeliv.tsv.gz")

metaG_presence_lessfiltered_simple_NULL <- simple_cooccur_parallel(in_df = metaG_presence_lessfiltered,
                                                                   ncores = 64,
                                                                   tmp_dir = "/mfs/gdouglas/tmp/metaG_presence_lessfiltered_simple_cooccur_tmp/",
                                                                   output_gzip_file = "/mfs/gdouglas/projects/ocean_mags/networks/filtersplit/metaG_simple_cooccur.lessfiltered.tsv.gz")
