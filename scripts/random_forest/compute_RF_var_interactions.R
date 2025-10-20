rm(list = ls(all.names = TRUE))

library(future)
library(iml)
library(ranger)
plan(multisession, workers = 8)

rf_output <- readRDS("/mfs/gdouglas/projects/ocean_hgt_zenodo/hgt_prev_analyses/random_forest_output/rf_output.rds")

null_based_df <- read.table("/mfs/gdouglas/projects/ocean_hgt_zenodo/hgt_prev_analyses/random_forest_output/null_varImp.tsv.gz",
                            sep = "\t", header = TRUE, row.names = 1)

spearman_info <- read.table("/mfs/gdouglas/projects/ocean_hgt_zenodo/hgt_prev_analyses/random_forest_output/var_spearman_vs_hgt_prop.tsv.gz",
                            sep = "\t", header = TRUE, row.names = 1)

breakdown <- read.table("/mfs/gdouglas/projects/ocean_hgt_zenodo/hgt_prev_analyses/Tara_lower_fraction_env_matched.tsv.gz",
                        sep = "\t", header = TRUE, stringsAsFactors = FALSE, row.names = 1)

cleaned_meta <- read.table(file = "/mfs/gdouglas/projects/ocean_hgt_zenodo/mapfiles/Tara_PANGEA_env_data.tsv.gz",
                           header=TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = 1)

cleaned_meta <- cleaned_meta[, -which(colnames(cleaned_meta) %in% c("Sample_ID_BioSamples_accession_number", "Sample_ID_ENA_sample_accession_number"))]

cleaned_meta <- cleaned_meta[rownames(breakdown), ]

cols_to_remove <- c("sample_w_diff_fraction",
                    "Fraction_lower_µm_used_on_board_to_prepare_samp",
                    "Fraction_upper_µm_used_on_board_to_prepare_samp",
                    "Sample_method_short_label_describing_the_ta",
                    "Station_TARA_event_datetime_station",
                    "Method_Device",
                    "Event",
                    "Date_Time",
                    "Sample_ID_TARA_barcode_registered_at",
                    "bac660_1_m",
                    "bacp_1_m")

cleaned_meta <- cleaned_meta[, -which(colnames(cleaned_meta) %in% cols_to_remove)]

cleaned_meta <- cleaned_meta[which(rowSums(is.na(cleaned_meta)) == 0), ]

variable_label_map <- list(
  "prop_hgt" = "prop_hgt",
  "depth_nominal" = "Depth",
  "Latitude" = "Latitude",
  "Longitude" = "Longitude",
  "Temperature" = "Temperature",
  "Salinity" = "Salinity",
  "Potential_density" = "Potential_density",
  "Oxygen_calibrated" = "Oxygen",
  "Nitrate" = "Nitrate",
  "Chlorophyll_a_NPQ_and_water_calib" = "Chlorophyll_a",
  "bbp470_1_m" = "bbp470",
  "fCDOM_ppb_QSE" = "fCDOM",
  "PAR_percent_8day_estimate" = "PAR")

colnames(cleaned_meta) <- sapply(colnames(cleaned_meta),
                                 function(x) variable_label_map[[x]])


null_based_df_sorted <- null_based_df[order(null_based_df$obs_permute_importance), ]
significant_vars <- rev(rownames(null_based_df_sorted)[null_based_df_sorted$null_p_BH < 0.05])

rf_input <- data.frame(prop_hgt = breakdown[rownames(cleaned_meta), 'prop_hgt'],
                       cleaned_meta)

predictor <- Predictor$new(rf_output, data = rf_input, y = 'prop_hgt')

pairwise_interactions <- data.frame(var1 = character(),
                                    var2 = character(),
                                    interaction_strength = numeric(),
                                    stringsAsFactors = FALSE)

for(i in 1:(length(significant_vars)-1)) {
  var1 <- significant_vars[i]
  var_interact <- Interaction$new(predictor, feature = var1)
  interact_results <- var_interact$results

  for(j in (i + 1):length(significant_vars)) {
    var2 <- significant_vars[j]

    interaction_strength <- interact_results[grep(var2, interact_results$.feature), '.interaction']

    if (nrow(interaction_strength) != 1) { stop('Error - no interaction found for this pair: ', var1, var2) }

    pairwise_interactions <- rbind(pairwise_interactions,
                                   data.frame(var1 = var1,
                                              var2 = var2,
                                              interaction_strength = as.numeric(interaction_strength)))
  }
}

write.table(pairwise_interactions,
            gzfile("/mfs/gdouglas/projects/ocean_hgt_zenodo/hgt_prev_analyses/random_forest_output/RF_sig_variable_pairwise_H_statistic.tsv.gz"),
            sep = "\t", row.names = FALSE, quote = FALSE)
