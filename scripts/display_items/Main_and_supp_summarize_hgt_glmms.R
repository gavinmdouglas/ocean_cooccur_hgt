rm(list = ls(all.names = TRUE))

library(ggplot2)
library(cowplot)

glmm_coefficients_raw <- list()

hgt_approaches <- c('clusterbased', 'blast', 'rangerdtl')
cooccur_approaches <- c('hyperg', 'simple', 'propr')

rds_filetypes <- list.files(path="/mfs/gdouglas/projects/ocean_mags/glmm_working/clusterbased/hyperg", pattern=".rds")

for (hgt in hgt_approaches) {
  for (cooccur in cooccur_approaches) {
    for (rds_file in rds_filetypes) {
      combo_id <- paste(hgt, cooccur, rds_file, sep = '_')
      glmm_coefficients_raw[[combo_id]] <- data.frame(readRDS(paste('/mfs/gdouglas/projects/ocean_mags/glmm_working/', hgt, cooccur, rds_file, sep = '/'))$coefficients$cond)
      glmm_coefficients_raw[[combo_id]]$BH <- p.adjust(glmm_coefficients_raw[[combo_id]]$Pr...z.., 'BH')
      glmm_coefficients_raw[[combo_id]] <- glmm_coefficients_raw[[combo_id]][which(glmm_coefficients_raw[[combo_id]]$BH < 0.05), ]
      glmm_coefficients_raw[[combo_id]]$cooccur_approach <- cooccur
      glmm_coefficients_raw[[combo_id]]$hgt_approach <- hgt
      glmm_coefficients_raw[[combo_id]]$glmm_type <- gsub('.rds$', '', rds_file)
      glmm_coefficients_raw[[combo_id]]$variable <- rownames(glmm_coefficients_raw[[combo_id]])
    }
  }
}

glmm_coefficients_all <- do.call(rbind, glmm_coefficients_raw)
rownames(glmm_coefficients_all) <- NULL



variable_map <- list()
variable_map[["(Intercept)"]] <- 'Intercept'
variable_map[["cooccur"]] <- 'Co-occurrence'
variable_map[["tip_dist_orderedNorm"]] <- 'Inter-tip distance'
variable_map[["depth_orderedNorm"]] <- 'Depth (diff.)'
variable_map[["latitude_orderedNorm"]] <- 'Latitude (diff.)'
variable_map[["longitude_orderedNorm"]] <- 'Longitude (diff.)'
variable_map[["temperature_orderedNorm"]] <- 'Temperature (diff.)'
variable_map[["oxygen_orderedNorm"]] <- 'Oxygen (diff.)'
variable_map[["salinity_orderedNorm"]] <- 'Salinity (diff.)'

var_order <- c('Co-occurrence', 'Inter-tip distance', 'Depth (diff.)', 'Latitude (diff.)', 'Longitude (diff.)', 'Temperature (diff.)', 'Oxygen (diff.)', 'Salinity (diff.)')

glmm_coefficients_all$clean_var <- NA
for (i in 1:nrow(glmm_coefficients_all)) {
  glmm_coefficients_all[i, 'clean_var'] <- variable_map[[glmm_coefficients_all$variable[i]]]
}

glmm_coefficients <- glmm_coefficients_all[which(glmm_coefficients_all$glmm_type == 'hgt_by_cooccur_tip_and_env'), ]

clusterbased_hyperg_glmm_coefficients <- glmm_coefficients[which(glmm_coefficients$hgt_approach == 'clusterbased' & glmm_coefficients$cooccur_approach == 'hyperg'), ]

clusterbased_hyperg_glmm_coefficients <- clusterbased_hyperg_glmm_coefficients[which(clusterbased_hyperg_glmm_coefficients$clean_var != 'Intercept'), ]

clusterbased_hyperg_glmm_coefficients$clean_var <- factor(clusterbased_hyperg_glmm_coefficients$clean_var,
                                                          levels = rev(var_order))

clusterbased_hyperg_glmm_coefficients$var_type <- c('Binary', rep('Continuous (ordered quantile normalized)', nrow(clusterbased_hyperg_glmm_coefficients) - 1))

clusterbased_hyperg_glmm_coefficients_barplot <- ggplot(data = clusterbased_hyperg_glmm_coefficients,
                                                        aes(x = Estimate, y = clean_var, fill = var_type)) +
  geom_bar(stat="identity") +
  scale_y_discrete(drop = FALSE) +

  theme_bw() +
  ylab("Variable") +
  xlab('Coefficient (log-odds)') +
  scale_fill_manual(name='Variable type', values=c('deepskyblue3', 'red3')) +
  geom_vline(xintercept = 0, linetype="dotted",
             color = "black") +
  geom_errorbar(aes(xmin = Estimate - Std..Error * 2,
                    xmax = Estimate + Std..Error * 2),
                width = 0.2, color = "black") +
  theme(legend.position = "inside",
        legend.position.inside = c(0.7, 0.3),
        legend.background = element_rect(fill = "white", color = "gray80"),
        legend.margin = margin(6, 6, 6, 6)) +
  ggtitle('Cluster-based HGT, HyperG co-occur (intercept: -18.02)')

ggsave(plot = clusterbased_hyperg_glmm_coefficients_barplot,
       filename = "/mfs/gdouglas/scripts/ocean_mag_hgt/display_items/Main_clusterbased_hyperg_glmm_coef.pdf",
       device = "pdf", width = 6, height = 4, units = "in", dpi=600)


# Same, but for "less-filtered" data.
glmm_coefficients_less_filtered <- glmm_coefficients_all[which(glmm_coefficients_all$glmm_type == 'hgt_by_cooccur_tip_and_env.lessfiltered'), ]

clusterbased_hyperg_glmm_coefficients <- glmm_coefficients_less_filtered[which(glmm_coefficients_less_filtered$hgt_approach == 'clusterbased' & glmm_coefficients_less_filtered$cooccur_approach == 'hyperg'), ]

intercept_tmp <- as.character(round(clusterbased_hyperg_glmm_coefficients[which(clusterbased_hyperg_glmm_coefficients$clean_var == 'Intercept'), 'Estimate'], 2))
clusterbased_hyperg_glmm_coefficients <- clusterbased_hyperg_glmm_coefficients[which(clusterbased_hyperg_glmm_coefficients$clean_var != 'Intercept'), ]

clusterbased_hyperg_glmm_coefficients$clean_var <- factor(clusterbased_hyperg_glmm_coefficients$clean_var,
                                                          levels = rev(var_order))

clusterbased_hyperg_glmm_coefficients$var_type <- c('Binary', rep('Continuous (ordered quantile)', nrow(clusterbased_hyperg_glmm_coefficients) - 1))

less_filtered_clusterbased_hyperg_glmm_coefficients_barplot <- ggplot(data = clusterbased_hyperg_glmm_coefficients,
                                                        aes(x = Estimate, y = clean_var, fill = var_type)) +
  geom_bar(stat="identity") +
  scale_y_discrete(drop = FALSE) +
  theme_bw() +
  ylab("Variable") +
  xlab('Coefficient (log-odds)') +
  scale_fill_manual(name='Variable type', values=c('deepskyblue3', 'red3')) +
  geom_vline(xintercept = 0, linetype="dotted",
             color = "black") +
  geom_errorbar(aes(xmin = Estimate - Std..Error * 2,
                    xmax = Estimate + Std..Error * 2),
                width = 0.2, color = "black") +
  theme(legend.position = "inside",
        legend.position.inside = c(0.75, 0.2),
        legend.background = element_rect(fill = "white", color = "gray80"),
        legend.text = element_text(size = 8),
       legend.title = element_text(size = 9)) +
  xlim(c(-3.5, 3.5)) +
  ggtitle(paste('Less-filtered (intercept: ', intercept_tmp, ')', sep = ''))


# Same, but for "freeliv" data.
glmm_coefficients_freeliv <- glmm_coefficients_all[which(glmm_coefficients_all$glmm_type == 'hgt_by_cooccur_tip_and_env.freeliv'), ]

clusterbased_hyperg_glmm_coefficients <- glmm_coefficients_freeliv[which(glmm_coefficients_freeliv$hgt_approach == 'clusterbased' & glmm_coefficients_freeliv$cooccur_approach == 'hyperg'), ]

intercept_tmp <- as.character(round(clusterbased_hyperg_glmm_coefficients[which(clusterbased_hyperg_glmm_coefficients$clean_var == 'Intercept'), 'Estimate'], 2))
clusterbased_hyperg_glmm_coefficients <- clusterbased_hyperg_glmm_coefficients[which(clusterbased_hyperg_glmm_coefficients$clean_var != 'Intercept'), ]

clusterbased_hyperg_glmm_coefficients$clean_var <- factor(clusterbased_hyperg_glmm_coefficients$clean_var,
                                                          levels = rev(var_order))

clusterbased_hyperg_glmm_coefficients$var_type <- c('Binary', rep('Continuous (ordered quantile)', nrow(clusterbased_hyperg_glmm_coefficients) - 1))

freeliv_clusterbased_hyperg_glmm_coefficients_barplot <- ggplot(data = clusterbased_hyperg_glmm_coefficients,
                                                                      aes(x = Estimate, y = clean_var, fill = var_type)) +
  geom_bar(stat="identity") +
  scale_y_discrete(drop = FALSE) +
  theme_bw() +
  ylab("Variable") +
  xlab('Coefficient (log-odds)') +
  scale_fill_manual(name='Variable type', values=c('deepskyblue3', 'red3')) +
  geom_vline(xintercept = 0, linetype="dotted",
             color = "black") +
  geom_errorbar(aes(xmin = Estimate - Std..Error * 2,
                    xmax = Estimate + Std..Error * 2),
                width = 0.2, color = "black") +
  theme(legend.position = "inside",
        legend.position.inside = c(0.75, 0.2),
        legend.background = element_rect(fill = "white", color = "gray80"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9)) +
  xlim(c(-3.5, 3.5)) +
  ggtitle(paste('Free-living (intercept: ' , intercept_tmp, ')', sep = ''))


filter_split_hyperg_clusterbased_glmm_combined <- plot_grid(freeliv_clusterbased_hyperg_glmm_coefficients_barplot,
                                                            less_filtered_clusterbased_hyperg_glmm_coefficients_barplot,
                                                            labels = c('a', 'b'),
                                                            nrow = 2)

ggsave(plot = filter_split_hyperg_clusterbased_glmm_combined,
       filename = "/mfs/gdouglas/scripts/ocean_mag_hgt/display_items/Main_filtersplit_clusterbased_hyperg_glmm_coef.pdf",
       device = "pdf", width = 6, height = 8, units = "in", dpi=600)


# Also make plots for all combinations and save to giant PDF, to explore.
pdf("/mfs/gdouglas/projects/ocean_mags/glmm_working/all_glmm_coef_plots.pdf")
for (hgt in hgt_approaches) {
  for (cooccur in cooccur_approaches) {
    for (rds_filetype in rds_filetypes) {
      rds_filetype <- gsub('.rds$', '', rds_filetype)
      tab_subset <- glmm_coefficients_all[which(glmm_coefficients_all$glmm_type == rds_filetype & glmm_coefficients_all$hgt_approach == hgt & glmm_coefficients_all$cooccur_approach == cooccur), ]

      intercept_tmp <- as.character(round(tab_subset[which(tab_subset$clean_var == 'Intercept'), 'Estimate'], 2))
      tab_subset <- tab_subset[which(tab_subset$clean_var != 'Intercept'), ]

      tab_subset$clean_var <- factor(tab_subset$clean_var, levels = rev(unique(tab_subset$clean_var)))

       plot_out <- ggplot(data = tab_subset, aes(x = Estimate, y = clean_var)) +
        geom_bar(stat="identity") +
        scale_y_discrete(drop = FALSE) +
        theme_bw() +
        ylab("Variable") +
        xlab('Coefficient (log-odds)') +
        geom_vline(xintercept = 0, linetype="dotted",
                   color = "black") +
        geom_errorbar(aes(xmin = Estimate - Std..Error * 2,
                          xmax = Estimate + Std..Error * 2),
                      width = 0.2, color = "black") +
        ggtitle(paste(hgt, cooccur, rds_filetype, '(intercept:' , paste(intercept_tmp, ')', sep = '')))

        print(plot_out)

    }
  }
}
dev.off()