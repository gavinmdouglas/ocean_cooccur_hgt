rm(list = ls(all.names = TRUE))

PRJEB9691_meta <- read.table('~/projects/ocean_mags/additional_protist_fraction/metadata/PRJEB9691_biosample_result.tsv',
                             header=TRUE, sep='\t', stringsAsFactors = FALSE, comment.char = '', quote='')
PRJEB9691_meta <- PRJEB9691_meta[which(PRJEB9691_meta$sampling_station != ''), ]
PRJEB9691_meta <- PRJEB9691_meta[which(PRJEB9691_meta$size_fraction_upper_threshold != '>'), ]
PRJEB9691_meta$size_fraction_upper_threshold <- as.numeric(PRJEB9691_meta$size_fraction_upper_threshold)

PRJEB9691_derep_raw <- list()
for (sampling_station in unique(PRJEB9691_meta$sampling_station)) {
  PRJEB9691_meta_station <- PRJEB9691_meta[which(PRJEB9691_meta$sampling_station == sampling_station), ]
  for (depth in unique(PRJEB9691_meta_station$depth)) {
    PRJEB9691_meta_station_depth <- PRJEB9691_meta_station[which(PRJEB9691_meta_station$depth == depth), ]
    PRJEB9691_meta_station_depth_focal <- PRJEB9691_meta_station_depth[which(PRJEB9691_meta_station_depth$size_fraction_lower_threshold == 3.0 & PRJEB9691_meta_station_depth$size_fraction_upper_threshold == 20), ]
    if (nrow(PRJEB9691_meta_station_depth_focal) == 0) { next }
    PRJEB9691_derep_raw[[paste(sampling_station, depth)]] <- PRJEB9691_meta_station_depth_focal[1, , drop = FALSE]
  }
}

PRJEB9691_derep <- do.call(rbind, PRJEB9691_derep_raw)

# PRJEB4352
PRJEB4352_meta <- read.table('~/projects/ocean_mags/additional_protist_fraction/metadata/PRJEB4352_biosample_result.tsv',
                             header=TRUE, sep='\t', stringsAsFactors = FALSE, comment.char = '', quote='')
PRJEB4352_meta <- PRJEB4352_meta[which(PRJEB4352_meta$sampling_station != ''), ]
PRJEB4352_meta <- PRJEB4352_meta[which(PRJEB4352_meta$size_fraction_upper_threshold != '>'), ]
PRJEB4352_meta$size_fraction_upper_threshold <- as.numeric(PRJEB4352_meta$size_fraction_upper_threshold)
PRJEB4352_meta <- PRJEB4352_meta[which(PRJEB4352_meta$size_fraction_lower_threshold != '>'), ]
PRJEB4352_meta <- PRJEB4352_meta[which(PRJEB4352_meta$size_fraction_lower_threshold != '0.8'), ]

PRJEB4352_derep_raw <- list()
for (sampling_station in unique(PRJEB4352_meta$sampling_station)) {
  PRJEB4352_meta_station <- PRJEB4352_meta[which(PRJEB4352_meta$sampling_station == sampling_station), ]
  for (depth in unique(PRJEB4352_meta_station$depth)) {
    PRJEB4352_meta_station_depth <- PRJEB4352_meta_station[which(PRJEB4352_meta_station$depth == depth), ]
    PRJEB4352_meta_station_depth_focal <- PRJEB4352_meta_station_depth[which(PRJEB4352_meta_station_depth$size_fraction_lower_threshold == 5 & PRJEB4352_meta_station_depth$size_fraction_upper_threshold == 20), ]
    if (nrow(PRJEB4352_meta_station_depth_focal) == 0) { next }
    PRJEB4352_derep_raw[[paste(sampling_station, depth)]] <- PRJEB4352_meta_station_depth_focal[1, , drop = FALSE]
  }
}
PRJEB4352_derep <- do.call(rbind, PRJEB4352_derep_raw)

PRJEB9691_derep$bioproject <- 'PRJEB9691'

PRJEB4352_derep$bioproject <- 'PRJEB4352'

# Add missing columns with NA values to each dataframe.
all_cols <- union(colnames(PRJEB9691_derep), colnames(PRJEB4352_derep))
for (col in all_cols) {
  if (! col %in% colnames(PRJEB9691_derep)) {
    PRJEB9691_derep[, col] <- NA
  }
  if (! col %in% colnames(PRJEB4352_derep)) {
    PRJEB4352_derep[, col] <- NA
  }
}

combined_filt <- rbind(PRJEB9691_derep[, all_cols], PRJEB4352_derep[, all_cols])

combined_filt[combined_filt == 99999] <- NA

write.table(x = combined_filt,
            file = '~/projects/ocean_mags/additional_protist_fraction/metadata/focal_derep_PRJEB9691_PRJEB4352_sample_meta.tsv',
            sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)


# Then to get run accessions to download.
PRJEB9691_run_info <- read.table('~/projects/ocean_mags/additional_protist_fraction/metadata/PRJEB9691_ENA_metadata.txt',
                                 sep = '\t', stringsAsFactors = FALSE, header=TRUE, comment.char = '', quote = '')

PRJEB9691_run_info_focal <- PRJEB9691_run_info[which(PRJEB9691_run_info$sample_accession %in% combined_filt$accession), ]

write.table(x = PRJEB9691_run_info_focal,
            file = '~/projects/ocean_mags/additional_protist_fraction/metadata/PRJEB9691_ENA_metadata_filt_focal.txt',
            sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)


PRJEB4352_run_info <- read.table('~/projects/ocean_mags/additional_protist_fraction/metadata/PRJEB4352_ENA_metadata.txt',
                                 sep = '\t', stringsAsFactors = FALSE, header=TRUE, comment.char = '', quote = '')

PRJEB4352_run_info_focal <- PRJEB4352_run_info[which(PRJEB4352_run_info$sample_accession %in% combined_filt$accession), ]

write.table(x = PRJEB4352_run_info_focal,
            file = '~/projects/ocean_mags/additional_protist_fraction/metadata/PRJEB4352_ENA_metadata_filt_focal.txt',
            sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)

# Also write out run accessions.
write.table(x = PRJEB9691_run_info_focal$run_accession,
            file = '~/projects/ocean_mags/additional_protist_fraction/metadata/run_accessions/PRJEB9691.txt',
            sep = '\t', col.names = FALSE, row.names = FALSE, quote = FALSE)

write.table(x = PRJEB4352_run_info_focal$run_accession,
            file = '~/projects/ocean_mags/additional_protist_fraction/metadata/run_accessions/PRJEB4352.txt',
            sep = '\t', col.names = FALSE, row.names = FALSE, quote = FALSE)

# Finally, get mapping of sample accessions to run accessions, to make it easier to concatenate technical replicates later.
sampleid_map_raw <- list()

for (sampleid in unique(PRJEB9691_run_info_focal$sample_accession)) {
  sample_subset <- PRJEB9691_run_info_focal[which(PRJEB9691_run_info_focal$sample_accession == sampleid), ]
  sampleid_map_raw[[sampleid]] <- data.frame(sample_accession = sampleid, run_ids = paste(sample_subset$run_accession, collapse = ';'))
}

for (sampleid in unique(PRJEB4352_run_info_focal$sample_accession)) {
  sample_subset <- PRJEB4352_run_info_focal[which(PRJEB4352_run_info_focal$sample_accession == sampleid), ]
  sampleid_map_raw[[sampleid]] <- data.frame(sample_accession = sampleid, run_ids = paste(sample_subset$run_accession, collapse = ';'))
}

sampleid_map <- do.call(rbind, sampleid_map_raw)

write.table(x = sampleid_map,
            file = '~/projects/ocean_mags/additional_protist_fraction/metadata/sample_accession_to_run_accession.tsv',
            sep = '\t', col.names = FALSE, row.names = FALSE, quote = FALSE)
