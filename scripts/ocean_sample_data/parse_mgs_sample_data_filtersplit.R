rm(list = ls(all.names = TRUE))

# Process all data for ocean metagenomics data related to environmental conditions, location, sampling time, etc.
# Outputs clean subset of columns of interest.
# Also, get mean values per sample.

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

# Then also get mean values for variable of interest per genome
# (based on the samples they are present across).
allsamples_presence <- read.table(file = "/mfs/gdouglas/projects/ocean_mags/networks/combined_tables/metaG_presence_allsamples.tsv.gz",
                                  header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)

rownames(info) <- info$sample_name

vars_of_interest <- c("lower_filter", "upper_filter", "depth", "latitude",
                      "longitude", "temperature", "oxygen", "salinity")

mean_tab <- data.frame(matrix(NA,
                              nrow = ncol(allsamples_presence),
                              ncol = length(vars_of_interest)))
colnames(mean_tab) <- vars_of_interest
rownames(mean_tab) <- colnames(allsamples_presence)

for (taxon in colnames(allsamples_presence)) {
  samples_w_taxon <- rownames(allsamples_presence)[which(allsamples_presence[, taxon] > 0)]
  if (length(samples_w_taxon) < 10) { stop("Error - should be at least 10 MGS samples!") }

  info_taxon_subset <- info[samples_w_taxon, ]

  for (var_of_interest in vars_of_interest) {
    vec <- info_taxon_subset[, var_of_interest]
    if (length(which(! is.na(vec))) > 0) {
      mean_tab[taxon, var_of_interest] <- mean(vec, na.rm = TRUE)
    }
  }
}

mean_gzfile_connection <-  gzfile("/mfs/gdouglas/projects/ocean_mags/networks/combined_tables/allsamples_present_metadata_mean_by_sample.tsv.gz", "w")
write.table(x = mean_tab,
            file = mean_gzfile_connection,
            col.names = NA, row.names = TRUE, sep = "\t", quote = FALSE)
close(mean_gzfile_connection)

median_tab <- data.frame(matrix(NA,
                              nrow = ncol(allsamples_presence),
                              ncol = length(vars_of_interest)))
colnames(median_tab) <- vars_of_interest
rownames(median_tab) <- colnames(allsamples_presence)

for (taxon in colnames(allsamples_presence)) {
  samples_w_taxon <- rownames(allsamples_presence)[which(allsamples_presence[, taxon] > 0)]
  if (length(samples_w_taxon) < 10) { stop("Error - should be at least 10 MGS samples!") }

  info_taxon_subset <- info[samples_w_taxon, ]

  for (var_of_interest in vars_of_interest) {
    vec <- info_taxon_subset[, var_of_interest]
    if (length(which(! is.na(vec))) > 0) {
      median_tab[taxon, var_of_interest] <- median(vec, na.rm = TRUE)
    }
  }
}

median_gzfile_connection <-  gzfile("/mfs/gdouglas/projects/ocean_mags/networks/combined_tables/allsamples_present_metadata_median_by_sample.tsv.gz", "w")
write.table(x = median_tab,
            file = median_gzfile_connection,
            col.names = NA, row.names = TRUE, sep = "\t", quote = FALSE)
close(median_gzfile_connection)
