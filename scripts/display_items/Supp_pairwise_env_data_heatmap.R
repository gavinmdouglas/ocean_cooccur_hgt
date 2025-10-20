rm(list = ls(all.names = TRUE))

library(ggplot2)
library(reshape2)
library(ComplexHeatmap)

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

cleaned_meta_scaled <- scale(cleaned_meta)
colnames(cleaned_meta_scaled) <- gsub('_', ' ', colnames(cleaned_meta_scaled))
cor_matrix <- cor(cleaned_meta, method = "spearman", use = "complete.obs")

# Create distance matrix based on absolute correlations
dist_matrix <- as.dist(1 - abs(cor_matrix))

env_heatmap <- Heatmap(t(cleaned_meta_scaled),
                      row_names_side = "right",
                      show_column_names = FALSE,
                      show_column_dend = FALSE,
                      row_dend_width = unit(3, "cm"),
                      clustering_distance_rows = dist_matrix,
                      clustering_distance_columns = "euclidean",
                      name = "Scaled\nvalue",
                      clustering_method_rows = "average",
                      clustering_method_columns = "average")

pdf("/mfs/gdouglas/scripts/ocean_cooccur_hgt/display_items/Supp_env_variable_heatmap.pdf", width = 9, height = 5)
env_heatmap
dev.off()
