rm(list = ls(all.names = TRUE))

# Parse Tara ocean metadata per sample, to subsample to SRRs and environmental factors of interest.

# Read in biosample IDs to keep:
PRJEB1787_info <- read.table('/mfs/gdouglas/projects/ocean_mags/metadata/PRJEB1787_metadata.csv', sep=',', stringsAsFactors = FALSE, header=TRUE)
PRJEB9740_info <- read.table('/mfs/gdouglas/projects/ocean_mags/metadata/PRJEB9740_metadata.csv', sep=',', stringsAsFactors = FALSE, header=TRUE)

protist_fraction <- read.table('/mfs/gdouglas/projects/ocean_mags/additional_protist_fraction/metadata/focal_derep_PRJEB9691_PRJEB4352_sample_meta.tsv',
                               sep = '\t', header=TRUE, stringsAsFactors = FALSE)

biosamples <- c(protist_fraction$accession, PRJEB1787_info$BioSample, PRJEB9740_info$BioSample)
biosamples <- unique(biosamples)

depth_env <- read.table('/mfs/gdouglas/projects/ocean_mags/metadata/Tara_env_data/dataset875582/datasets/cleaned/TARA_ENV_DEPTH_SENSORS.txt',
                        header=TRUE, sep = '\t', comment.char = "")

depth_env <- depth_env[which(depth_env$Sample.ID..BioSamples.accession.number...... %in% biosamples), ]

depth_env_initial <- depth_env[, 1:24]

depth_env_medians <- depth_env[, grep("median", colnames(depth_env))]

# Manually select some columns and rename clearly, based on descriptions of them in the raw file from PANGAEA
# (confusingly they have the same column names, but their order matters).

# Note that I would have liked to keep the in situ PAR estimates, but these were mainly missing (NA).
colnames(depth_env_medians)[which(colnames(depth_env_medians) == "PAR..mol.quanta.m..2.day...median.value..50th.percentile.....2")] <- 'PAR_absolute_8day_estimate'
colnames(depth_env_medians)[which(colnames(depth_env_medians) == "PAR......median.value..50th.percentile.....2")] <- 'PAR_percent_8day_estimate'

colnames(depth_env_medians)[which(colnames(depth_env_medians) == "Chl.a..mg.m..3...median.value..50th.percentile.....2")] <- 'Chlorophyll_a_NPQ_and_water_calib'

colnames(depth_env_medians)[which(colnames(depth_env_medians) == "beta470..m.sr...median.value..50th.percentile.....2")] <- 'beta470_dark_calibrated'

colnames(depth_env_medians)[which(colnames(depth_env_medians) == "O2..µmol.kg...median.value..50th.percentile.....1")] <- 'Oxygen_calibrated'

colnames(depth_env_medians)[which(colnames(depth_env_medians) == "X.NO3....µmol.l...median.value..50th.percentile....")] <- 'Nitrate'


median_col_to_rm <- c("PAR..mol.quanta.m..2.day...median.value..50th.percentile....", "PAR......median.value..50th.percentile....",
                      "PAR..mol.quanta.m..2.day...median.value..50th.percentile.....1", "PAR..mol.quanta.m..2.day...median.value..50th.percentile.....3",
                      "PAR......median.value..50th.percentile.....1", "PAR......median.value..50th.percentile.....3",
                      "Chl.a..mg.m..3...median.value..50th.percentile....", "Chl.a..mg.m..3...median.value..50th.percentile.....1",
                      "beta470..m.sr...median.value..50th.percentile....", "beta470..m.sr...median.value..50th.percentile.....1",
                      "O2..µmol.kg...median.value..50th.percentile....")

depth_env_medians <- depth_env_medians[, -which(colnames(depth_env_medians) %in% median_col_to_rm)]

# Create base names by removing the trailing numbers
colnames(depth_env_medians) <- gsub("\\.\\.+\\d+$", "", colnames(depth_env_medians))

colnames(depth_env_medians) <- gsub("\\.median\\.value.*", "", colnames(depth_env_medians))

colnames(depth_env_medians)[which(colnames(depth_env_medians) == "Temp...C..")] <- 'Temperature'
colnames(depth_env_medians)[which(colnames(depth_env_medians) == "Cond..mS.cm..")] <- 'Conductivity'
colnames(depth_env_medians)[which(colnames(depth_env_medians) == "Sal.")] <- 'Salinity'
colnames(depth_env_medians)[which(colnames(depth_env_medians) == "Tpot...C..")] <- 'Potential_temperature'
colnames(depth_env_medians)[which(colnames(depth_env_medians) == "Sigma.theta..kg.m..3..")] <- "Potential_density"

colnames(depth_env_initial)[which(colnames(depth_env_initial) == "Depth.nominal..from.which.this.sample.was.co....")] <- "depth_nominal"

depth_env_initial$depth_nominal[depth_env_initial$depth_nominal == "5-160"] <- mean(c(5, 160))
depth_env_initial$depth_nominal[depth_env_initial$depth_nominal == "10-100"] <- mean(c(10, 100))
depth_env_initial$depth_nominal <- as.numeric(depth_env_initial$depth_nominal)

initial_col_to_rm <- c("Campaign",
                       "Basis",
                       "Distance..km...interval.between.the.location....",
                       "Duration..interval.between.the.date.tim....",
                       "File.name..of.sensor.profiles.that.meet.....",
                       "NOBS......that.meet.the.criteria..10.m.....",
                       "Sample.label..TARA_event.datetime_station._....",
                       "Samp.mat..TARA_station._environmental.f....",
                       "Env.feature...abbreviation...full.name..EN....",
                       "Depth.top..m...from.which.this.sample.was.co....",
                       "Depth.bot..m...from.which.this.sample.was.co....")

depth_env_initial <- depth_env_initial[, -which(colnames(depth_env_initial) %in% initial_col_to_rm)]

focal_depth_env <- cbind(depth_env_initial, depth_env_medians)

colnames(focal_depth_env) <- gsub("\\.", "_", colnames(focal_depth_env))
colnames(focal_depth_env) <- gsub("___", "_", colnames(focal_depth_env))
colnames(focal_depth_env) <- gsub("__", "_", colnames(focal_depth_env))
colnames(focal_depth_env) <- gsub("_*$", "", colnames(focal_depth_env))

# Now check for very similar variables, to remove
pairwise_spearman_df_cont <- function(df) {
  df <- df[sapply(df, is.numeric)]
  pairs <- t(combn(colnames(df), 2))
  result <- data.frame(Var1 = pairs[, 1],
                       Var2 = pairs[, 2],
                       cor = NA,
                       pvalue = NA)

  for(i in 1:nrow(pairs)) {
    test <- cor.test(df[[pairs[i, 1]]], df[[pairs[i, 2]]], method = "spearman", exact=FALSE)
    result$cor[i] <- test$estimate
    result$pvalue[i] <- test$p.value
  }

  return(result)
}

focal_depth_env_cor <- pairwise_spearman_df_cont(focal_depth_env)

focal_depth_env_cor_0.95 <- focal_depth_env_cor[abs(focal_depth_env_cor$cor) >= 0.95, ]

# Kept only bbp470_1_m and dropped other optic properties variables that were highly correlated.
# Similarly, remove Conductivity and Potential_temperature (redundant with Temperature)
# Also PAR_absolute_8day_estimate redundant with PAR_percent_8day_estimate
redundant_to_rm <- c("beta470_dark_calibrated", "bb470_1_m", "Conductivity", "Potential_temperature", "PAR_absolute_8day_estimate")

focal_depth_env <- focal_depth_env[, -which(colnames(focal_depth_env) %in% redundant_to_rm)]


# Finally identify samples with identical environmental features, which should not be treated independently.
# These correspond to cases where the same sample was size fractionated using different cut-offs.
identify_group_num <- function(df) {
  group_ids <- rep(NA, nrow(df))
  group_counter <- 1

  # Handle rows that are not duplicated
  non_dup_rows <- which(!duplicated(df) & !duplicated(df, fromLast = TRUE))
  for (row_idx in non_dup_rows) {
    group_ids[row_idx] <- paste0("group", group_counter)
    group_counter <- group_counter + 1
  }

  # Handle duplicated rows
  rows_to_consider <- which(duplicated(df) | duplicated(df, fromLast = TRUE))
  if (length(rows_to_consider) == 0) { return(group_ids) }

  df_dup <- df[rows_to_consider, ]
  df_unique <- df_dup[-which(duplicated(df_dup)), ]

  for (i in 1:nrow(df_unique)) {
    group_id <- paste0("group", group_counter)
    rows_to_rm <- integer()
    for (j in rows_to_consider) {
      if (all(df[j, ] == df_unique[i, ], na.rm = TRUE)) {
        group_ids[j] <- group_id
        rows_to_rm <- c(rows_to_rm, j)
      }
    }
    rows_to_consider <- setdiff(rows_to_consider, rows_to_rm)
    group_counter <- group_counter + 1
  }
  return(group_ids)
}

# First for rows that are 100% identical, excluding sample IDs.
# Sanity check, as there should be no technical replicates.
focal_depth_env_sample_rm <- focal_depth_env[, 4:ncol(focal_depth_env)]
focal_depth_env_sample_rm <- focal_depth_env_sample_rm[, -which(colnames(focal_depth_env_sample_rm) == "Sample_method_short_label_describing_the_ta")]
all_env_identical_groups <- identify_group_num(focal_depth_env_sample_rm)
length(unique(all_env_identical_groups)) == nrow(focal_depth_env_sample_rm)

# Then based on all features except fraction (i.e. filter) cut-offs
focal_depth_env_sample_fraction_rm <- focal_depth_env_sample_rm[, -grep("^Fraction", colnames(focal_depth_env_sample_rm))]
fraction_diff_groups <- identify_group_num(focal_depth_env_sample_fraction_rm)
length(unique(fraction_diff_groups)) == nrow(focal_depth_env_sample_rm)

# Also run a check that none of the samples are identical simply when depth is ignored.
focal_depth_env_sample_depth_rm <- focal_depth_env_sample_rm[, -grep("^depth_nominal", colnames(focal_depth_env_sample_rm))]
depth_diff_groups <- identify_group_num(focal_depth_env_sample_depth_rm)
length(unique(depth_diff_groups)) == nrow(focal_depth_env_sample_depth_rm)

initial_first3col <- colnames(focal_depth_env)[1:3]
initial_remaining_col <- colnames(focal_depth_env)[4:ncol(focal_depth_env)]

focal_depth_env$sample_w_diff_fraction <- fraction_diff_groups


# And then to get samples with at least 1 genome found in our dataset, and add this sample name as first column.
# I foolishly mixed different sample ID types with the protist-fraction and OceanDNA samples, so this takes an extra step
# to figure out the right sample ID to use.
prev <- read.table('/mfs/gdouglas/projects/ocean_mags/networks/combined_tables/metaG_presence_tara_w_protist_fraction.tsv.gz',
                   header=TRUE, sep = '\t', stringsAsFactors = FALSE, row.names = 1)

focal_depth_env$sample_match <- rep(NA, nrow(focal_depth_env))
focal_depth_env$sample_match[which(focal_depth_env$Sample_ID_BioSamples_accession_number %in% rownames(prev))] <- focal_depth_env$Sample_ID_BioSamples_accession_number[which(focal_depth_env$Sample_ID_BioSamples_accession_number %in% rownames(prev))]
focal_depth_env$sample_match[which(focal_depth_env$Sample_ID_ENA_sample_accession_number %in% rownames(prev))] <- focal_depth_env$Sample_ID_ENA_sample_accession_number[which(focal_depth_env$Sample_ID_ENA_sample_accession_number %in% rownames(prev))]
focal_depth_env <- focal_depth_env[which(! is.na(focal_depth_env$sample_match)), ]


new_col_order <- c("sample_match", initial_first3col, "sample_w_diff_fraction", initial_remaining_col)
focal_depth_env <- focal_depth_env[, new_col_order]

write.table(x=focal_depth_env, file=gzfile("/mfs/gdouglas/projects/ocean_hgt_zenodo/mapfiles/Tara_PANGEA_env_data.tsv.gz"),
            sep="\t", col.names = TRUE, row.names = FALSE, quote = TRUE)
