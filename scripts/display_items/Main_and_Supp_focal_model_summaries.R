rm(list = ls(all.names = TRUE))

library(ComplexHeatmap)

model_info <- readRDS(file = '/mfs/gdouglas/projects/ocean_mags/glmm_working/combined_summary.rds')
names(model_info) <- gsub('//', '/', names(model_info))
vif_info <- readRDS('/mfs/gdouglas/projects/ocean_mags/glmm_working/all_vif.rds')
names(vif_info) <- gsub('//', '/', names(vif_info))

coef_vars <- rownames(model_info$`/mfs/gdouglas/projects/ocean_mags/glmm_working/tara_genomes/blast_allsamples/hyperg/hgt_by_cooccur_tip_and_env_w_filtergroup.glm.rds`$coef)

details <- data.frame(matrix(NA, nrow = length(model_info), ncol = 38))
colnames(details) <- c("genome_subset", "tool_and_subset", "approach", "model", "max_vif", coef_vars, paste(coef_vars, "CI", sep = '_'), paste(coef_vars, "P", sep = '_'))
rownames(details) <- gsub('.rds', '', names(model_info))

for (i in 1:length(model_info)) {
  details[i, "genome_subset"] <- model_info[[i]]$genomes
  details[i, "tool_and_subset"] <- model_info[[i]]$tool_and_subset
  details[i, "approach"] <- model_info[[i]]$approach
  details[i, "model"] <- sub('.rds', '', basename(names(model_info)[i]))

  i_coef <- numeric()
  i_coef_CI <- character()
  i_coef_p <- numeric()

  coef_tab <- model_info[[i]]$coef

  if (length(grep('.glmm.rds', basename(names(model_info)[i]))) > 0) {
    coef_tab <- coef_tab$cond
  }

  for (coef in coef_vars) {
    if (coef %in% rownames(coef_tab)) {
      i_coef <- c(i_coef, coef_tab[coef, 'Estimate'])
      lower_ci <- sprintf("%.2f", coef_tab[coef, 'Estimate'] - 1.96 * coef_tab[coef, 'Std. Error'])
      upper_ci <- sprintf("%.2f", coef_tab[coef, 'Estimate'] + 1.96 * coef_tab[coef, 'Std. Error'])
      i_coef_CI <- c(i_coef_CI, paste(lower_ci, upper_ci, sep = ','))
      i_coef_p <- c(i_coef_p, coef_tab[coef, 'Pr(>|z|)'])
    } else {
      i_coef <- c(i_coef, NA)
      i_coef_CI <- c(i_coef_CI, NA)
      i_coef_p <- c(i_coef_p, NA)
    }
  }

  details[i, coef_vars] <- i_coef
  details[i, paste(coef_vars, "CI", sep = '_')] <- i_coef_CI
  details[i, paste(coef_vars, "P", sep = '_')] <- i_coef_p

  exp_vif_file <-  sub('.rds', '.VIF.tsv', names(model_info)[i])
  if (exp_vif_file %in% names(vif_info)) {
    details[i, 'max_vif'] <- max(vif_info[[exp_vif_file]]$vif$VIF)
  }
}

# write.table(x = details, file = '/mfs/gdouglas/projects/ocean_hgt_zenodo/hgt_cooccur_enrich/model_summaries.tsv',
#             quote=FALSE, row.names = FALSE, col.names = TRUE, sep = '\t')


# Make cleaned-up table/heatmap with subset of models that are most informative.

focal_model_i <- c(which(details$genome_subset == 'tara_genomes' & details$tool_and_subset == "clusterbased_allsamples" & details$model == 'hgt_by_cooccur_tip_and_env_w_filtergroup.glm'),
                   which(details$genome_subset == 'tara_genomes' & details$tool_and_subset == "clusterbased_allsamples" & details$model == 'hgt_by_cooccur_tip_and_env_w_taxon_w_filtergroup.glmm'),

                   which(details$genome_subset == 'tara_genomes' & details$tool_and_subset == "clusterbased_geotraces_samples_hyperg" & details$model == 'hgt_by_cooccur_tip_and_env.glm'),
                   which(details$genome_subset == 'tara_genomes' & details$tool_and_subset == "clusterbased_tara_samples_hyperg" & details$model == 'hgt_by_cooccur_tip_and_env.glm'),

                   which(details$genome_subset == 'tara_genomes' & details$tool_and_subset == "blast_allsamples" & details$model == 'hgt_by_cooccur_tip_and_env_w_filtergroup.glm'),
                   which(details$genome_subset == 'tara_genomes' & details$tool_and_subset == "blast_allsamples" & details$model == 'hgt_by_cooccur_tip_and_env_w_taxon_w_filtergroup.glmm'),
                   which(details$genome_subset == 'tara_genomes' & details$tool_and_subset == "rangerdtl_allsamples" & details$model == 'hgt_by_cooccur_tip_and_env_w_filtergroup.glm'),
                   which(details$genome_subset == 'tara_genomes' & details$tool_and_subset == "rangerdtl_allsamples" & details$model == 'hgt_by_cooccur_tip_and_env_w_taxon_w_filtergroup.glmm'),
                   which(details$genome_subset == 'progenomes' & details$tool_and_subset == "clusterbased_hyperg_allsamples" & details$model == 'hgt_by_cooccur_tip_and_env_w_filtergroup.glm'),
                   which(details$genome_subset == 'progenomes' & details$tool_and_subset == "clusterbased_hyperg_geotraces_samples" & details$model == 'hgt_by_cooccur_tip_and_env.glm'),
                   which(details$genome_subset == 'progenomes' & details$tool_and_subset == "clusterbased_hyperg_tara_samples" & details$model == 'hgt_by_cooccur_tip_and_env.glm'))

details_focal <- details[focal_model_i, ]
details_focal_to_plot <- details_focal

clean_var <- character()
for (variable in coef_vars) {
  clean_colname <- paste(variable, 'clean', sep = '_')
  ci_colname <- paste(variable, 'CI', sep = '_')
  p_colname <- paste(variable, 'P', sep = '_')
  details_focal[, clean_colname] <- NA
  clean_var <- c(clean_var, clean_colname)
  for (j in 1:nrow(details_focal)) {
    if (is.na(details_focal[j, variable])) {
      details_focal[j, clean_colname] <- '-'
    } else if (details_focal[j, p_colname] >= 0.05) {
      details_focal[j, clean_colname] <- 'n.s.'
    } else {
      #details_focal[j, clean_colname] <- paste(sprintf("%.2f", details_focal[j, variable]), '\n(', details_focal[j, ci_colname], ')', sep = '')
      details_focal[j, clean_colname] <- sprintf("%.2f", details_focal[j, variable])
    }

  }

  details_focal_to_plot[which(details_focal_to_plot[, p_colname] >= 0.05), variable] <- NA
  details_focal_to_plot[which(details_focal_to_plot[, variable] > 0), variable] <- 'Positive'
  details_focal_to_plot[which(details_focal_to_plot[, variable] < 0), variable] <- 'Negative'

}

details_focal_to_plot$max_vif <- NA

details_focal_to_plot$genome_subset <- gsub('tara_genomes', 'Focal genomes', details_focal_to_plot$genome_subset)
details_focal_to_plot$genome_subset <- gsub('progenomes', 'proGenomes', details_focal_to_plot$genome_subset)

details_focal_to_plot$tool_and_subset[grep('clusterbased', details_focal_to_plot$tool_and_subset)] <- 'Cluster-based HGT'
details_focal_to_plot$tool_and_subset[grep('blast', details_focal_to_plot$tool_and_subset)] <- 'BLAST-based HGT'
details_focal_to_plot$tool_and_subset[grep('rangerdtl', details_focal_to_plot$tool_and_subset)] <- 'RANGER-DTL-based HGT'

details_focal_to_plot$approach[grep('hyperg', details_focal_to_plot$approach)] <- 'HyperG co-occur.'
details_focal_to_plot$approach[grep('propr', details_focal_to_plot$approach)] <- 'propr co-abun.'
details_focal_to_plot$approach[grep('simple', details_focal_to_plot$approach)] <- 'Simple co-occur.'



model_rownames <- sapply(1:nrow(details_focal_to_plot), function(i) { paste(details_focal_to_plot$genome_subset[i], details_focal_to_plot$tool_and_subset[i], details_focal_to_plot$approach[i], sep = ' | ') })

for (i in which(details_focal_to_plot$model == "hgt_by_cooccur_tip_and_env_w_taxon_w_filtergroup.glmm")) {
  model_rownames[i] <- paste(model_rownames[i], "(GLMM)")
}

model_rownames[grep('tara_genomes/clusterbased_geotraces_samples_hyperg/hgt_by_cooccur_tip_and_env.glm', rownames(details_focal_to_plot))] <-  "Focal genomes | Cluster-based HGT | HyperG co-occur. (GEOTRACES samples only)"
model_rownames[grep('tara_genomes/clusterbased_tara_samples_hyperg/hgt_by_cooccur_tip_and_env.glm', rownames(details_focal_to_plot))] <-  "Focal genomes | Cluster-based HGT | HyperG co-occur. (Tara samples only)"

model_rownames[grep('progenomes/clusterbased_hyperg_geotraces_samples/hgt_by_cooccur_tip_and_env.glm', rownames(details_focal_to_plot))] <-  "proGenomes | Cluster-based HGT | HyperG co-occur. (GEOTRACES samples only)"
model_rownames[grep('progenomes/clusterbased_hyperg_tara_samples/hgt_by_cooccur_tip_and_env.glm', rownames(details_focal_to_plot))] <-  "proGenomes | Cluster-based HGT | HyperG co-occur. (Tara samples only)"

model_rownames[grep('progenomes/clusterbased_hyperg_allsamples/hgt_by_cooccur_tip_and_env_w_filtergroup.glm', rownames(details_focal_to_plot))] <-  "proGenomes | Cluster-based HGT | HyperG co-occur."

model_breakdown <- details_focal_to_plot[, c(1:3)]
colnames(model_breakdown) <- c('Genome DB', 'HGT method', 'Co-occur method')
model_breakdown$`Co-occur method`[grep('.rds$', model_breakdown$`Co-occur method`)] <- 'HyperG'
model_breakdown$`Co-occur method` <- gsub(' co-occur.', '', model_breakdown$`Co-occur method`)
model_breakdown$`Co-occur method` <- gsub(' co-abun.', '', model_breakdown$`Co-occur method`)

details_focal_text <- details_focal[, c(clean_var, "max_vif")]

details_focal_text$max_vif <- sprintf("%.2f", details_focal_text$max_vif)

details_focal_to_plot <- details_focal_to_plot[, c(coef_vars, "max_vif")]

supp_row_i <- c(grep('.glmm$', rownames(details_focal_to_plot)), grep("_geotraces_", rownames(details_focal_to_plot)), grep("_tara_", rownames(details_focal_to_plot)))

details_focal_text_supp <- details_focal_text[supp_row_i, ]
details_focal_to_plot_supp <- details_focal_to_plot[supp_row_i, ]
model_breakdown_supp <- model_breakdown[supp_row_i, ]
model_rownames_supp <- model_rownames[supp_row_i]

details_focal_text <- details_focal_text[-supp_row_i, ]
details_focal_to_plot <- details_focal_to_plot[-supp_row_i, ]
model_breakdown <- model_breakdown[-supp_row_i, ]
model_rownames <- model_rownames[-supp_row_i]

rownames(details_focal_text) <- model_rownames
rownames(details_focal_to_plot) <- model_rownames

genome_annot <- HeatmapAnnotation(df = model_breakdown,
                                     col = list(`Genome DB` = c("Focal genomes" = "grey80", "proGenomes" = "black"),
                                         `HGT method` = c("Cluster-based HGT" = "dodgerblue1", "BLAST-based HGT" = "lightskyblue1", "RANGER-DTL-based HGT" = "blue2"),
                                         `Co-occur method` = c("HyperG" = "green3", propr = "palegreen1", "Simple" = "green4")))

details_focal_to_plot <- t(details_focal_to_plot)
details_focal_text <- t(details_focal_text)

col_to_plot <- c('indianred1', 'skyblue1')
names(col_to_plot) <- c('Positive', 'Negative')
main_heatmap <- ComplexHeatmap::Heatmap(matrix = details_focal_to_plot,
                                        col=col_to_plot,
                                        na_col = 'grey95',
                                        top_annotation = genome_annot,
                                        cluster_rows = FALSE,
                                        cluster_columns = FALSE,
                                        cell_fun = function(j, i, x, y, width, height, fill) {
                                          grid.text(details_focal_text[i, j], x, y, gp = gpar(fontsize = 10))
                                        },
                                        name = 'Association',
                                        row_names_centered = FALSE,
                                        column_title = 'Model',
                                        row_split = c(rep('', 11), ' '),
                                        show_column_names = FALSE,
                                        row_labels = c('Intercept', 'Co-occur.', 'Tip dist.', "Filter group ('Free-living')", "Filter group ('Less-filtered')", 'Depth', 'Latitude', 'Longitude', 'Temperature', 'Oxygen', 'Salinity', 'Max. VIF')
)

pdf("~/scripts/ocean_mag_hgt/display_items/Main_all_model_summary.pdf", width = 8, height = 5)
draw(main_heatmap, padding = unit(c(2, 2, 2, 2), "mm"), heatmap_legend_side = "top", annotation_legend_side = "top",  merge_legend=TRUE)
dev.off()



rownames(details_focal_text_supp) <- model_rownames_supp
rownames(details_focal_to_plot_supp) <- model_rownames_supp

genome_annot_supp <- HeatmapAnnotation(df = model_breakdown_supp,
                                  col = list(`Genome DB` = c("Focal genomes" = "grey80", "proGenomes" = "black"),
                                             `HGT method` = c("Cluster-based HGT" = "dodgerblue1", "BLAST-based HGT" = "lightskyblue1", "RANGER-DTL-based HGT" = "blue2"),
                                             `Co-occur method` = c("HyperG" = "green3", propr = "palegreen1", "Simple" = "green4")))

details_focal_to_plot_supp <- t(details_focal_to_plot_supp)
details_focal_text_supp <- t(details_focal_text_supp)

col_split_factor <- factor(c(rep('GLMM (on duplicated table) with taxon random effect', 9), rep('GEOTRACES\nsamples only', 2), rep('Tara\nsamples only', 2)),
                           levels = c('GLMM (on duplicated table) with taxon random effect', 'GEOTRACES\nsamples only', 'Tara\nsamples only'))


supp_heatmap <- ComplexHeatmap::Heatmap(matrix = details_focal_to_plot_supp,
                                        col=col_to_plot,
                                        na_col = 'grey95',
                                        top_annotation = genome_annot_supp,
                                        cluster_rows = FALSE,
                                        cluster_columns = FALSE,
                                        cell_fun = function(j, i, x, y, width, height, fill) {
                                          grid.text(details_focal_text_supp[i, j], x, y, gp = gpar(fontsize = 10))
                                        },
                                        name = 'Association',
                                        row_names_centered = FALSE,
                                        #column_title = 'Model',
                                        row_split = c(rep('', 11), ' '),
                                        show_column_names = FALSE,
                                        row_labels = c('Intercept', 'Co-occur.', 'Tip dist.', "Filter group ('Free-living')", "Filter group ('Less-filtered')", 'Depth', 'Latitude', 'Longitude', 'Temperature', 'Oxygen', 'Salinity', 'Max. VIF'),
                                        column_split =col_split_factor
)

pdf("~/scripts/ocean_mag_hgt/display_items/Supp_model_summary.pdf", width = 10, height = 5)
draw(supp_heatmap, padding = unit(c(2, 2, 2, 2), "mm"), heatmap_legend_side = "top", annotation_legend_side = "top",  merge_legend=TRUE)
dev.off()

