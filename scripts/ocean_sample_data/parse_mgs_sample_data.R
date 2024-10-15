rm(list = ls(all.names = TRUE))

# Process all data for ocean metagenomics data related to environmental conditions, location, sampling time, etc.
# Outputs clean subset of columns of interest.
# Also, get mean values per sample.

source('~/scripts/ocean_mag_hgt/scripts/ocean_sample_data/compute_mean_and_median_samplevals_function.R')

info <- read.table("/mfs/gdouglas/projects/ocean_mags/metadata/OceanDNA_supp_metadata/Supp_File_S1_water_samples.tsv",
                       header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Re-code one sample that for some reason is listed as two IDs (just use the one in the metadata file downloaded from SRA).
info[which(info$sample_name == "ERS492821_ERS492814"), "sample_name"] <- "ERS492814"

colnames(info)[which(colnames(info) == "longigute")] <- "longitude"

# Variables of interest:
# depth
# latitude and longitude
# temperature
# oxygen
# salinity

info <- info[, c("sample_name", "bioproject", "instrument", "collection_date",
                 "depth", "latitude", "longitude", "lower_filter",
                 "upper_filter", "temperature", "oxygen", "salinity"), ]

# Replace all hyphens with NA
info[info == "-"] <- NA

# Convert several columns wrongly encoded as chars.
info$lower_filter <- as.numeric(info$lower_filter)
info$upper_filter <- as.numeric(info$upper_filter)
info$temperature <- as.numeric(info$temperature)
info$oxygen <- as.numeric(info$oxygen)
info$salinity <- as.numeric(info$salinity)

# Write out this table.
write.table(x = info,
            file = "/mfs/gdouglas/projects/ocean_mags/metadata/OceanDNA_supp_metadata/Supp_File_S1_water_samples_clean.tsv",
            col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

rownames(info) <- info$sample_name

# Then also get mean values for variable of interest per genome
# (based on the samples they are present across).
allsamples_presence <- read.table(file = "/mfs/gdouglas/projects/ocean_mags/networks/combined_tables/metaG_presence_allsamples.tsv.gz",
                                  header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)

compute_mean_median_sample_vals(info_tab=info,
                                presence_tab=allsamples_presence,
                                outprefix='/mfs/gdouglas/projects/ocean_mags/networks/combined_tables/allsamples_present_metadata')

# Also get mean and median values for filter-split samples.
metaG_presence_freeliv <- read.table(file = "/mfs/gdouglas/projects/ocean_mags/networks/combined_tables/metaG_presence_freeliv.tsv.gz",
                                     header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)

metaG_presence_lessfiltered <- read.table(file = "/mfs/gdouglas/projects/ocean_mags/networks/combined_tables/metaG_presence_lessfiltered.tsv.gz",
                                          header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)

compute_mean_median_sample_vals(info_tab=info,
                                presence_tab=metaG_presence_freeliv,
                                outprefix='/mfs/gdouglas/projects/ocean_mags/networks/combined_tables/freeliv_present_metadata')

compute_mean_median_sample_vals(info_tab=info,
                                presence_tab=metaG_presence_lessfiltered,
                                outprefix='/mfs/gdouglas/projects/ocean_mags/networks/combined_tables/lessfiltered_present_metadata')
