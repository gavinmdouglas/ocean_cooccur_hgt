rm(list = ls(all.names = TRUE))

library(ggplot2)
library(reshape2)
library(iml)
library(ranger)
library(cowplot)

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

sig_env <- cleaned_meta[, significant_vars]

sig_env$prop_hgt <- breakdown[rownames(cleaned_meta), "prop_hgt"]

sig_env_long <- reshape2::melt(data=sig_env, id.vars="prop_hgt")
sig_env_long$variable <- factor(sig_env_long$variable, levels=significant_vars)

spearman_info$rho_label <- paste("Spearman =", round(spearman_info$rho, 3))

sig_varImp_spearman_info <- spearman_info[significant_vars, ]
sig_varImp_spearman_info$variable <- rownames(sig_varImp_spearman_info)
sig_varImp_spearman_info$variable <- factor(sig_varImp_spearman_info$variable, levels=sig_varImp_spearman_info$variable)

sig_varImp_spearman_info$hjust <- ifelse(rownames(sig_varImp_spearman_info) %in% c("Nitrate", "Depth", "fCDOM"), -0.1, 1.1)
sig_varImp_spearman_info$x_pos <- ifelse(rownames(sig_varImp_spearman_info) %in% c("Nitrate", "Depth", "fCDOM"), -Inf, Inf)

RF_sig_var_scatterplots <- ggplot(sig_env_long, aes(x = value, y = prop_hgt)) +
  geom_point(alpha = 0.6) +
  geom_text(data = sig_varImp_spearman_info, aes(x = x_pos, y = Inf, label = rho_label, hjust = hjust),
            vjust = 1.5, size = 3, inherit.aes = FALSE, colour="grey70") +
  facet_wrap(~ variable, scales = "free_x", labeller = as_labeller(function(x) gsub("_", " ", x))) +
  labs(x = "Environmental variable value", y = "Horizontal gene transfer prevalence") +
  theme_bw()


# Also get Individual Conditional Expectation plots for each significant variable.
rf_input <- data.frame(prop_hgt = breakdown[rownames(cleaned_meta), 'prop_hgt'],
                       cleaned_meta)

predictor <- Predictor$new(rf_output, data = rf_input, y = 'prop_hgt')

ice_results_raw <- list()

for (sig_var in significant_vars) {
  ice_out <- FeatureEffect$new(predictor, feature = sig_var, method = "ice")
  ice_out$results$Variable <- sig_var
  colnames(ice_out$results)[1] <- "x_val"
  ice_results_raw[[sig_var]] <- ice_out$results
}

ice_combined <- do.call(rbind, ice_results_raw)

rf_input_sig_var <- reshape2::melt(rf_input[, significant_vars],
                                   variable.name = "Variable",
                                   value.name = "value")
rf_input_sig_var$Variable <- as.character(rf_input_sig_var$Variable)

ice_combined$Variable <- factor(ice_combined$Variable, levels = significant_vars)

rf_input_sig_var$Variable <- factor(rf_input_sig_var$Variable, levels = significant_vars)

ice_linegraphs <- ggplot(data = ice_combined, aes(x = x_val, y = .value)) +
  geom_line(aes(group = .id), alpha = 0.3, color = "grey") +
  stat_summary(aes(group = 1), fun = mean, geom = "line",
               color = "blue", linewidth = 1) +
  geom_rug(data = rf_input_sig_var, aes(x = value),
           inherit.aes = FALSE, alpha = 0.5) +  # Rug plot
  facet_wrap(~ Variable, scales = "free_x", labeller = labeller(Variable = function(x) gsub("_", " ", x))) +
  labs(x = "Environmental variable value", y = "Predicted horizontal gene transfer prevalence") +
  theme_bw()



# Plot pairwise interactions for significant variables.
pairwise_interactions <- read.table("/mfs/gdouglas/projects/ocean_hgt_zenodo/hgt_prev_analyses/random_forest_output/RF_sig_variable_pairwise_H_statistic.tsv.gz",
                                    sep = "\t", stringsAsFactors = FALSE, header=TRUE)
pairwise_interactions$var1 <- gsub("_", " ", pairwise_interactions$var1)

interaction_matrix <- matrix(NA, nrow = length(significant_vars), ncol = length(significant_vars))
rownames(interaction_matrix) <-  gsub("_", " ", significant_vars)
colnames(interaction_matrix) <- gsub("_", " ", significant_vars)

for(i in 1:nrow(pairwise_interactions)) {
  var1 <- pairwise_interactions$var1[i]
  var2 <- pairwise_interactions$var2[i]
  value <- pairwise_interactions$interaction_strength[i]

  interaction_matrix[var1, var2] <- value
  interaction_matrix[var2, var1] <- value
}

heatmap_data <- melt(interaction_matrix)
names(heatmap_data) <- c("var1", "var2", "interaction")

heatmap_data <- heatmap_data[as.numeric(heatmap_data$var1) > as.numeric(heatmap_data$var2), ]

heatmap_data$var1 <- factor(heatmap_data$var1,
                            levels = gsub("_", " ", significant_vars))

heatmap_data$var2 <- factor(heatmap_data$var2,
                            levels = gsub("_", " ", significant_vars))

interaction_plot <- ggplot(heatmap_data, aes(x = var2, y = var1, fill = interaction)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red", limits=c(0, 0.125)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "", y = "", fill = "H-statistic")

combined_indepth_RF_eval_lower_row <- plot_grid(NULL, interaction_plot, NULL, rel_widths = c(0.55, 1, 0.55), ncol=3, labels=c('', 'b', ''))
grey_line <- ggplot() +
  theme_void() +
  theme(plot.background = element_rect(fill = "grey", color = "grey"))

combined_indepth_RF_eval_plots <- plot_grid(ice_linegraphs,
                                            grey_line,
                                            combined_indepth_RF_eval_lower_row,
                                            nrow=3,
                                            rel_heights = c(1, 0.01, 0.5),
                                            labels=c('a', '', ''))

ggsave(filename = "/mfs/gdouglas/scripts/ocean_cooccur_hgt/display_items/Supp_sig_RF_var_scatterplots.pdf",
       plot = RF_sig_var_scatterplots, height=6, width=8, units = "in", device="pdf", dpi=300)

ggsave(filename = "/mfs/gdouglas/scripts/ocean_cooccur_hgt/display_items/Supp_sig_RF_var_indepth_info.pdf",
       plot = combined_indepth_RF_eval_plots, height=8, width=8, units = "in", device="pdf", dpi=300)
