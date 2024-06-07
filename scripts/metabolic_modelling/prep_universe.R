rm(list = ls(all.names = TRUE))

# Process all data for ocean metagenomics data related to environmental conditions, location, sampling time, etc.
# This script simply subsets the table to samples overlapping with the "all-sample" CoverM output, and
# outputs a clean table for downstream use.

TableS1 <- read.table("/mfs/gdouglas/projects/ocean_mags/metadata/OceanDNA_supp_metadata/Supp_File_S1.tsv",
                      header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Re-code one sample that for some reason is listed as two IDs (just use the one in the metadata file downloaded from SRA).
TableS1[which(TableS1$sample_name == "ERS492821_ERS492814"), "sample_name"] <- "ERS492814"

allsamples_presence <- read.table("/mfs/gdouglas/projects/ocean_mags/coverm/combined_tables/metaG_presence_allsamples.tsv.gz",
                                  header = TRUE, sep = "\t", row.names = 1)

info <- TableS1[which(TableS1$sample_name %in% rownames(allsamples_presence)), ]

if (length(setdiff(rownames(allsamples_presence), info$sample_name)) > 0 | length(setdiff(info$sample_name, rownames(allsamples_presence))) > 0) {
  stop("Not all samples matched") 
}

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
            file = "/mfs/gdouglas/projects/ocean_mags/coverm/combined_tables/allsamples_present_metadata.tsv",
            col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

# Then also get mean values for variable of interest per genome
# (based on the samples they are present across).
rownames(info) <- info$sample_name

vars_of_interest <- c("depth", "latitude", "longitude", "temperature", "oxygen", "salinity")

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

write.table(x = mean_tab,
            file = "/mfs/gdouglas/projects/ocean_mags/coverm/combined_tables/allsamples_present_metadata_mean_by_sample.tsv",
            col.names = NA, row.names = TRUE, sep = "\t", quote = FALSE)
