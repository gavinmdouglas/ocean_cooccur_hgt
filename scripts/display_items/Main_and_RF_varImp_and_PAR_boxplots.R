rm(list = ls(all.names = TRUE))

library(ggplot2)
library(reshape2)
library(cowplot)
library(ggbeeswarm)

null_based_df <- read.table("/mfs/gdouglas/projects/ocean_hgt_zenodo/hgt_prev_analyses/random_forest_output/null_varImp.tsv.gz",
                            sep = "\t", header = TRUE, row.names = 1)

rf_output <- readRDS("/mfs/gdouglas/projects/ocean_hgt_zenodo/hgt_prev_analyses/random_forest_output/rf_output.rds")

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
  "fCDOM_ppb_QSE" = "fCDOM",
  "PAR_percent_8day_estimate" = "PAR",
  "bbp470_1_m" = "bbp470")

colnames(cleaned_meta) <- sapply(colnames(cleaned_meta),
                                 function(x) variable_label_map[[x]])

varimp_plot_df <- data.frame(variable = names(rf_output$variable.importance),
                             importance = rf_output$variable.importance,
                             correlation = spearman_info[names(rf_output$variable.importance), "rho"],
                             p_adj = null_based_df[names(rf_output$variable.importance), "null_p_BH"],
                             significant = null_based_df[names(rf_output$variable.importance), "null_p_BH"] < 0.05)

varimp_plot_df$color <- ifelse(!varimp_plot_df$significant, "grey",
                               ifelse(varimp_plot_df$correlation > 0, "red", "blue"))

varimp_plot_df$association <- factor(
  ifelse(! varimp_plot_df$significant, "Not significant",
         ifelse(varimp_plot_df$correlation > 0, "Positive", "Negative")),
  levels = c("Positive", "Negative", "Not significant"))

varimp_plot_df <- varimp_plot_df[order(varimp_plot_df$importance, decreasing = FALSE), ]
varimp_plot_df$variable <- factor(varimp_plot_df$variable, levels = varimp_plot_df$variable)

varimp_plot <- ggplot(varimp_plot_df, aes(x = variable, y = importance, fill = association)) +
                geom_col() +
                scale_fill_manual(values = c("Positive" = "red",
                                             "Negative" = "blue",
                                             "Not significant" = "grey"),
                                  name = "") +
                coord_flip() +
                labs(x = "Variable",
                     y = "Permutation-based variable importance\n(for predicting horizontal gene transfer prevalence)") +
                theme_bw() +
                theme(axis.text.y = element_text(size = 10),
                      legend.position = "bottom") +
                scale_x_discrete(labels = function(x) gsub("_", " ", x))

# Then get panel of PAR vs hgt_prop split by extreme PAR values.
low_PAR_samples <- rownames(cleaned_meta)[which(cleaned_meta$PAR < 20)]
high_PAR_samples <- rownames(cleaned_meta)[which(cleaned_meta$PAR >= 20)]

PAR_vs_hgt_df <- data.frame(prop_hgt=breakdown[rownames(cleaned_meta), 'prop_hgt'],
                            PAR_grouping=NA)
rownames(PAR_vs_hgt_df) <- rownames(cleaned_meta)

PAR_vs_hgt_df[low_PAR_samples, 'PAR_grouping'] <- 'Low PAR'
PAR_vs_hgt_df[high_PAR_samples, 'PAR_grouping'] <- 'High PAR'

PAR_vs_hgt_df$PAR_grouping <- factor(PAR_vs_hgt_df$PAR_grouping, levels=c('Low PAR', 'High PAR'))

mean_diff <- mean(PAR_vs_hgt_df$prop_hgt[PAR_vs_hgt_df$PAR_grouping == "Low PAR"]) -
             mean(PAR_vs_hgt_df$prop_hgt[PAR_vs_hgt_df$PAR_grouping == "High PAR"])
wilcox_p <- wilcox.test(prop_hgt ~ PAR_grouping, data = PAR_vs_hgt_df)$p.value

PAR_vs_hgt_plot <- ggplot(data=PAR_vs_hgt_df, aes(x = PAR_grouping, y = prop_hgt)) +
  geom_quasirandom(alpha = 0.6, width = 0.3) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA, width = 0.5) +
  labs(x = "Sample grouping", y = "Horizontal gene transfer prevalence") +
  annotate("text", x = 1.5, y = 0.23,
           label = paste0("Mean difference: ", round(mean_diff, 4), "\n Wilcoxon p < 0.001"),
           hjust = -0.1, vjust = 1.1, size = 3.5) +
  theme_bw()


main_combined <- plot_grid(varimp_plot, PAR_vs_hgt_plot, labels=c('a', 'b'), rel_widths = c(1, 1))

ggsave(filename = "/mfs/gdouglas/scripts/ocean_cooccur_hgt/display_items/Main_Figure6.pdf",
       plot = main_combined, height=4, width=9.5, units = "in", device="pdf", dpi=600)
