rm(list = ls(all.names = TRUE))

library(ranger)
library(ggplot2)
library(parallel)

# Make sure that RF variable importance results don't change much if one of a pair of highly correlated variables are left out
# (the worry is that the importance is being split across multiple correlated variables, so they seem less predictive than they are in reality)

run_rf_analysis <- function(exclude_vars = NULL, num_trees = 10000, n_permutations = 1000, mc_cores = 200) {

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

  colnames(cleaned_meta) <- sapply(colnames(cleaned_meta), function(x) variable_label_map[[x]])


  if (!is.null(exclude_vars)) {
    cleaned_meta <- cleaned_meta[, -which(colnames(cleaned_meta) %in% exclude_vars)]
  }

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

  rownames(spearman_info) <- spearman_info$var

  rf_input <- data.frame(prop_hgt = breakdown[rownames(cleaned_meta), 'prop_hgt'], cleaned_meta)

  set.seed(97026)
  rf_ranger <- ranger(prop_hgt ~ ., data = rf_input, importance = "permutation", num.trees = num_trees)

  run_null_permutation <- function(i) {
    shuffled_data <- rf_input
    shuffled_data$prop_hgt <- sample(shuffled_data$prop_hgt)
    rf_null <- ranger(shuffled_data$prop_hgt ~ ., data = shuffled_data, num.trees = num_trees, importance = "permutation")

    return(list(
      importance = rf_null$variable.importance[names(rf_ranger$variable.importance)],
      r_squared = rf_null$r.squared
    ))
  }

  results <- mclapply(1:n_permutations, run_null_permutation, mc.cores = mc_cores)

  null_importance <- data.frame(matrix(NA, nrow = n_permutations, ncol = length(rf_ranger$variable.importance)))
  colnames(null_importance) <- names(rf_ranger$variable.importance)
  null_rsquared <- numeric(n_permutations)

  for (i in 1:n_permutations) {
    null_importance[i, ] <- results[[i]]$importance
    null_rsquared[i] <- results[[i]]$r_squared
  }

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

  varimp_plot_df <- data.frame(variable = names(rf_ranger$variable.importance),
                               importance = rf_ranger$variable.importance,
                               correlation = spearman_info[names(rf_ranger$variable.importance), "rho"],
                               p_adj = null_based_df[names(rf_ranger$variable.importance), "null_p_BH"],
                               significant = null_based_df[names(rf_ranger$variable.importance), "null_p_BH"] < 0.05)

  varimp_plot_df$association <- factor(
    ifelse(! varimp_plot_df$significant, "Not significant",
           ifelse(varimp_plot_df$correlation > 0, "Positive", "Negative")),
    levels = c("Positive", "Negative", "Not significant"))

  varimp_plot_df <- varimp_plot_df[order(varimp_plot_df$importance, decreasing = FALSE), ]
  varimp_plot_df$variable <- factor(varimp_plot_df$variable, levels = varimp_plot_df$variable)

  varimp_plot <- ggplot(varimp_plot_df, aes(x = variable, y = importance, fill = association)) +
    geom_col() +
    scale_fill_manual(values = c("Positive" = "red", "Negative" = "blue", "Not significant" = "grey"), name = "") +
    coord_flip() +
    labs(x = "Variable", y = "Permutation-based variable importance\n(for predicting horizontal gene transfer prevalence)") +
    theme_bw() +
    theme(axis.text.y = element_text(size = 10), legend.position = "bottom") +
    scale_x_discrete(labels = function(x) gsub("_", " ", x))

  return(list(rf_model = rf_ranger, null_results = null_based_df, spearman_results = spearman_info, plot = varimp_plot))
}

# Run analyses with different variable exclusions
no_depth_test <- run_rf_analysis(exclude_vars = "Depth")
no_PAR_test <- run_rf_analysis(exclude_vars = "PAR")
no_depth_PAR_test <- run_rf_analysis(exclude_vars = c("Depth", "PAR"))

no_chl_test <- run_rf_analysis(exclude_vars = "Chlorophyll_a")
no_bbp470_test <- run_rf_analysis(exclude_vars = "bbp470")
no_chlo_bbp470_test <- run_rf_analysis(exclude_vars = c("Chlorophyll_a", "bbp470"))

no_temperature_test <- run_rf_analysis(exclude_vars = "Temperature")
no_potential_density_test <- run_rf_analysis(exclude_vars = "Potential_density")
no_temperature_or_potential_density_test <- run_rf_analysis(exclude_vars = c("Temperature", "Potential_density"))

no_nitrate_test <- run_rf_analysis(exclude_vars = "Nitrate")
no_fCDOM_test <- run_rf_analysis(exclude_vars = "fCDOM")
no_nitrate_or_fCDOM_test <- run_rf_analysis(exclude_vars = c("Nitrate", "fCDOM"))


no_depth_test$plot
no_PAR_test$plot
no_depth_PAR_test$plot

no_chl_test$plot
no_bbp470_test$plot
no_chlo_bbp470_test$plot

no_temperature_test$plot
no_potential_density_test$plot
no_temperature_or_potential_density_test$plot

no_nitrate_test$plot
no_fCDOM_test$plot
no_nitrate_or_fCDOM_test$plot
