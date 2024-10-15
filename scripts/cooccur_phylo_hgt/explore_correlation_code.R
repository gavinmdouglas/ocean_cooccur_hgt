approaches <- c("simple", "hyperG", "propr")

env_vars <- c("depth", "latitude", "longitude", "temperature", "oxygen", "salinity")

id_map = list()

for (env_var in env_vars) {
  id_map[[env_var]] <- paste(tools::toTitleCase(env_var), "(diff.)")
}

id_map[["tip_dist"]] <- "Tip dist."

id_map[["both_gene_count"]] <- "HGT (gene count)"
id_map[["ranger_hgt_tallies"]] <- "HGT (gene count)"

id_map[["cooccur_simple_cooccur"]] <- "Co-occur (Simple)"
id_map[["cooccur_ratio"]] <- "Co-occur (HyperG)"
id_map[["cooccur_asso"]] <- "Co-occur (propr)"

approach_to_colname <- list()
approach_to_colname[["simple"]] <- "cooccur_simple_cooccur"
approach_to_colname[["hyperG"]] <- "cooccur_ratio"
approach_to_colname[["propr"]] <- "cooccur_asso"

clean_ids <- function(in_vec) {
  sapply(in_vec, function(x) { id_map[[x]] })
}

comparison_type_cols <- c("#33a02c", "#1f78b4", "#e31a1c")

summary_stats <- function(x) {
  n <- length(x)
  mean_val <- mean(x, na.rm = TRUE)
  median_val <- median(x, na.rm = TRUE)
  sd_val <- sd(x, na.rm = TRUE)
  se_val <- sd_val / sqrt(n)
  wilcox_out <- wilcox.test(x, exact=FALSE)
  out <- c(n, mean_val, median_val, sd_val, se_val, wilcox_out$statistic, wilcox_out$p.value)
  names(out) <- c("n", "mean", "median", "sd", "se", "Wilcox_V", "Wilcox_P")
  return(out)
}

read_corr_info <- function(pairwise_file, partial_file) {
  pairwise_in <- read.table(pairwise_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  pairwise_in$sig <- "Non-sig."
  pairwise_in$sig[which(pairwise_in$p.unc < 0.05)] <- "P < 0.05"

  partial_in <- read.table(partial_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

  return(list(pairwise_tab=pairwise_in,
              partial_tab=partial_in))
}

prep_pairwise_tab <- function(infile) {
  in_pairwise <- read.table(infile, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  X_tallies <- table(in_pairwise$X)
  Y_tallies <- table(in_pairwise$Y)
  in_pairwise$X <- factor(in_pairwise$X, levels = names(X_tallies)[order(X_tallies, decreasing = FALSE)])
  in_pairwise$Y <- factor(in_pairwise$Y, levels = names(Y_tallies)[order(Y_tallies, decreasing = FALSE)])
  in_pairwise_mean <- aggregate(x = r ~ X + Y, FUN = mean, data = in_pairwise)
  in_pairwise_mean_wide <- reshape2::dcast(data = in_pairwise_mean, formula = X ~ Y, value.var = 'r')
  rownames(in_pairwise_mean_wide) <- in_pairwise_mean_wide[, 1]

  # Also get comparisons not significantly different from 0, so that these can be coloured as grey.
  unique_combos <- in_pairwise[, c("X", "Y")]
  unique_combos <- unique_combos[-which(duplicated(unique_combos)), ]
  nonsig_comparisons <- list()
  nonsig_count <- 1
  for (i in 1:nrow(unique_combos)) {
    x_subset <- unique_combos[i, "X"]
    y_subset <- unique_combos[i, "Y"]
    combo_subset <- in_pairwise[which(in_pairwise$X == x_subset & in_pairwise$Y == y_subset), ]
    wilox_combo_out <- wilcox.test(combo_subset$r, exact=FALSE)
    if (wilox_combo_out$p.value >= 0.05) {
      # print("nonsig")
      nonsig_comparisons[[nonsig_count]] <- c(y_subset, x_subset)
      nonsig_count <- nonsig_count + 1
    }
  }

  return(list(tab=in_pairwise_mean_wide[, -1],
              nonsig_comparisons=nonsig_comparisons))
}

read_in_pairwise <- function(subfolder, approaches, suffix, ref_approach='simple') {
  mean_pairwise = list()

  for (approach in approaches) {
    infile_blast_pairwise <- paste(subfolder, approach, suffix, sep = "")
    mean_pairwise[[approach]] <- prep_pairwise_tab(infile_blast_pairwise)
  }

  identity_sanity_checks <- c()
  for (approach1 in approaches) {
    for (approach2 in approaches) {
      if (approach1 == approach2) { next }
      identity_sanity_checks <- c(identity_sanity_checks, identical(mean_pairwise[[approach1]]$tab[-8, ],  mean_pairwise[[approach2]]$tab[-8, ]))
    }
  }
  if (length(which(! identity_sanity_checks)) > 0) { stop("Non-co-occurrence mismatches" )}

  # Also make it easy to combine the co-occurrence types into a single heatmap.
  mean_pairwise_combined <- mean_pairwise[[ref_approach]]$tab
  for (approach in approaches) {
    if (approach == ref_approach) { next }
    tmp_rows <- mean_pairwise[[approach]]$tab[approach_to_colname[[approach]], , drop = FALSE]
    mean_pairwise_combined <- rbind(mean_pairwise_combined, tmp_rows)

  }

  return(list(mean_pairwise=mean_pairwise,
              mean_pairwise_combined=mean_pairwise_combined))
}


create_markdown_pairwise_heatmaps <- function(mean_pairwise_in, in_mat, approaches, category) {

  cat("\n\n")

  cat("### ", category, " results", "\n\n")

  in_mat <- as.matrix(in_mat)

  in_mat_char <- in_mat
  for (i in 1:ncol(in_mat)) {
    in_mat_char[, i] <- sprintf("%.2f", in_mat[, i])
  }
  in_mat_char[in_mat_char == "NA"] <- ""
  col_fun = colorRamp2(c(-1, 0, 1), c("#6699CC", "white", "#CC3333"))

  for (approach in approaches) {
    if (length(mean_pairwise_in[[approach]]$nonsig_comparisons) > 0) {
      for (non_sig_element in mean_pairwise_in[[approach]]$nonsig_comparisons) {
        non_sig_element <- as.character(non_sig_element)
        if (non_sig_element[2] %in% colnames(mean_pairwise_in[[approach]]$tab) & non_sig_element[1] %in% rownames(mean_pairwise_in[[approach]]$tab)) {
          in_mat[non_sig_element[1], non_sig_element[2]] <- NA
        }

        if (non_sig_element[1] %in% colnames(mean_pairwise_in[[approach]]$tab) & non_sig_element[2] %in% rownames(mean_pairwise_in[[approach]]$tab)) {
          in_mat[non_sig_element[2], non_sig_element[1]] <- NA
        }

      }

    }
  }

  heatmap_out <- ComplexHeatmap::Heatmap(matrix = in_mat,
                                         col = col_fun,
                                         cluster_rows = FALSE,
                                         cluster_columns = FALSE,
                                         row_split = c(" ", " ", " ", " ", " ", "  ", "   ", "     ", "     ", "     "),
                                         row_labels = clean_ids(rownames(in_mat)),
                                         column_labels = clean_ids(colnames(in_mat)),
                                         na_col = "grey50",
                                         cell_fun = function(j, i, x, y, width, height, fill) {
                                           grid.text(in_mat_char[i, j], x, y, gp = gpar(fontsize = 10))
                                         },
                                         name = 'Spearman correlation')

  print(heatmap_out)

  cat("\n\n")

}

