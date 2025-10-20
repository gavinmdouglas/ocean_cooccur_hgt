rm(list = ls(all.names = TRUE))

library(ranger)
library(ggplot2)

breakdown <- read.table("/mfs/gdouglas/projects/ocean_hgt_zenodo/hgt_prev_analyses/Tara_lower_fraction_env_matched.tsv.gz",
                        sep = "\t", header = TRUE, stringsAsFactors = FALSE, row.names = 1)

cleaned_meta <- read.table(file = "/mfs/gdouglas/projects/ocean_hgt_zenodo/mapfiles/Tara_PANGEA_env_data.tsv.gz",
                           header=TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = 1)

cleaned_meta <- cleaned_meta[, -which(colnames(cleaned_meta) %in% c("Sample_ID_BioSamples_accession_number", "Sample_ID_ENA_sample_accession_number"))]

cleaned_meta <- cleaned_meta[rownames(breakdown), ]


# Remove metadata columns.
# Also, decided later to also remove the scatter-related variables based on a
# recommendation that they were too hard to interpret (and including them)
# didn't change the overall signal.

cols_to_remove <- c("sample_w_diff_fraction",
                    "Fraction_lower_µm_used_on_board_to_prepare_samp",
                    "Fraction_upper_µm_used_on_board_to_prepare_samp",
                    "Sample_method_short_label_describing_the_ta",
                    "Station_TARA_event_datetime_station",
                    "Method_Device",
                    "Event",
                    "Date_Time",
                    "Sample_ID_TARA_barcode_registered_at",
                    'bac660_1_m',
                    'bacp_1_m')

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
  "fCDOM_ppb_QSE" = "fCDOM",
  "PAR_percent_8day_estimate" = "PAR",
  "bbp470_1_m" = "bbp470")

colnames(cleaned_meta) <- sapply(colnames(cleaned_meta),
                             function(x) variable_label_map[[x]])

cleaned_meta <- cleaned_meta[, sort(colnames(cleaned_meta))]

spearman_cor <- numeric()
spearman_p <- numeric()
for (varname in colnames(cleaned_meta)) {
  cor_out <- cor.test(cleaned_meta[, varname], breakdown[rownames(cleaned_meta), 'prop_hgt'], method='spearman')
  spearman_cor <- c(spearman_cor, cor_out$estimate)
  spearman_p <- c(spearman_p, cor_out$p.value)
}
spearman_info <- data.frame(var = colnames(cleaned_meta),
                            rho = spearman_cor,
                            p = spearman_p)

rf_input <- data.frame(prop_hgt = breakdown[rownames(cleaned_meta), 'prop_hgt'],
                       cleaned_meta)

num_trees_set <- 10000

rf_ranger <- ranger(prop_hgt ~ .,
                    data = rf_input,
                    importance = "permutation",
                    num.trees = num_trees_set,
                    seed = 97026)

run_null_permutation <- function(i) {
  shuffled_data <- rf_input
  shuffled_data$prop_hgt <- sample(shuffled_data$prop_hgt)
  rf_null <- ranger(prop_hgt ~ .,
                    data = shuffled_data,
                    num.trees = num_trees_set,
                    importance = "permutation",
                    seed = i)

  return(list(
    importance = rf_null$variable.importance[names(rf_ranger$variable.importance)],
    r_squared = rf_null$r.squared
  ))
}

results <- parallel::mclapply(1:1000, run_null_permutation, mc.cores = 200)

null_importance <- data.frame(matrix(NA, nrow = 1000, ncol = length(rf_ranger$variable.importance)))
colnames(null_importance) <- names(rf_ranger$variable.importance)
null_rsquared <- numeric(1000)

for (i in 1:1000) {
  null_importance[i, ] <- results[[i]]$importance
  null_rsquared[i] <- results[[i]]$r_squared
}

# Compute P-value based Monte Carlo permutation test.
null_based_p <- numeric()
for (varname in names(rf_ranger$variable.importance)) {
  numerator <- length(which(null_importance[, varname] >= rf_ranger$variable.importance[varname])) + 1
  denominator <- nrow(null_importance) + 1
  null_based_p <- c(null_based_p, numerator / denominator)
}
names(null_based_p) <- names(rf_ranger$variable.importance)

null_based_bh <- p.adjust(null_based_p, 'BH')

null_based_df <- data.frame(variable=names(null_based_bh),
                            obs_permute_importance=rf_ranger$variable.importance[names(null_based_bh)],
                            mean_null_permute_importance=colMeans(null_importance)[names(null_based_bh)],
                            raw_null_p=null_based_p,
                            null_p_BH=null_based_bh)

write.table(null_based_df,
            gzfile("/mfs/gdouglas/projects/ocean_hgt_zenodo/hgt_prev_analyses/random_forest_output/null_varImp.tsv.gz"),
            sep = "\t", row.names = FALSE, quote = FALSE)

saveRDS(rf_ranger, "/mfs/gdouglas/projects/ocean_hgt_zenodo/hgt_prev_analyses/random_forest_output/rf_output.rds")

write.table(spearman_info,
            gzfile("/mfs/gdouglas/projects/ocean_hgt_zenodo/hgt_prev_analyses/random_forest_output/var_spearman_vs_hgt_prop.tsv.gz"),
            sep = "\t", row.names = FALSE, quote = FALSE)
