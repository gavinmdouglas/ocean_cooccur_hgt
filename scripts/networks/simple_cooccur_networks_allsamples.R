rm(list = ls(all.names = TRUE))

source("~/scripts/ocean_mag_hgt/scripts/networks/network_functions.R")

metaG_presence <- read.table(file = "/mfs/gdouglas/projects/ocean_mags/networks/combined_tables/metaG_presence_allsamples.tsv.gz",
                             header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)

metaG_simple_NULL <- simple_cooccur_parallel(in_df = metaG_presence,
                                             ncores = 64,
                                             tmp_dir = "/mfs/gdouglas/tmp/metaG_simple_cooccur_tmp/",
                                             output_gzip_file = "/mfs/gdouglas/projects/ocean_mags/networks/allsamples/metaG_simple_cooccur.allsamples.tsv.gz")
