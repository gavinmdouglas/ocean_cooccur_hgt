rm(list = ls(all.names = TRUE))

library(ggplot2)
library(cowplot)

model_info <- readRDS(file = '/mfs/gdouglas/projects/ocean_mags/glmm_working/combined_summary.rds')

focal_model <- model_info$`/mfs/gdouglas/projects/ocean_mags/glmm_working/tara_genomes//clusterbased_allsamples/hyperg/hgt_by_cooccur_tip_and_env_w_filtergroup.glm.rds`

vif_in <- read.table('/mfs/gdouglas/projects/ocean_mags/glmm_working/tara_genomes/clusterbased_allsamples/hyperg/hgt_by_cooccur_tip_and_env_w_filtergroup.glm.VIF.tsv',
                     header=TRUE, sep='\t', stringsAsFactors = FALSE)

clusterbased_hyperg_coef <- data.frame(focal_model$coef)
clusterbased_hyperg_coef$BH <- p.adjust(clusterbased_hyperg_coef$Pr...z.., 'BH')
clusterbased_hyperg_coef$variable <- rownames(clusterbased_hyperg_coef)

variable_map <- list()
variable_map[["(Intercept)"]] <- 'Intercept'
variable_map[["cooccur"]] <- 'Co-occurrence'
variable_map[["tip_dist_orderedNorm"]] <- 'Inter-tip distance'
variable_map[["filter_group_matchFree-living"]] <- '\"Free-living\" enriched'
variable_map[["filter_group_matchLess-filtered"]] <- '\"Less filtered\" enriched'
variable_map[["depth_orderedNorm"]] <- 'Depth (diff.)'
variable_map[["latitude_orderedNorm"]] <- 'Latitude (diff.)'
variable_map[["longitude_orderedNorm"]] <- 'Longitude (diff.)'
variable_map[["temperature_orderedNorm"]] <- 'Temperature (diff.)'
variable_map[["oxygen_orderedNorm"]] <- 'Oxygen (diff.)'
variable_map[["salinity_orderedNorm"]] <- 'Salinity (diff.)'

var_order <- c('Co-occurrence', '\"Free-living\" enriched', '\"Less filtered\" enriched', 'Inter-tip distance', 'Depth (diff.)', 'Latitude (diff.)', 'Longitude (diff.)', 'Temperature (diff.)', 'Oxygen (diff.)', 'Salinity (diff.)')

clusterbased_hyperg_coef$clean_var <- NA
for (i in 1:nrow(clusterbased_hyperg_coef)) {
  clusterbased_hyperg_coef[i, 'clean_var'] <- variable_map[[clusterbased_hyperg_coef$variable[i]]]
}

clusterbased_hyperg_coef <- clusterbased_hyperg_coef[which(clusterbased_hyperg_coef$clean_var != 'Intercept'), ]

clusterbased_hyperg_coef$clean_var <- factor(clusterbased_hyperg_coef$clean_var,
                                                          levels = rev(var_order))

rownames(clusterbased_hyperg_coef) <- clusterbased_hyperg_coef$clean_var
clusterbased_hyperg_coef <- clusterbased_hyperg_coef[var_order, ]

clusterbased_hyperg_coef$var_type <- c(rep('Binary', 3), rep('Continuous', nrow(clusterbased_hyperg_coef) - 3))

clusterbased_hyperg_coef_barplot <- ggplot(data = clusterbased_hyperg_coef,
                                                        aes(x = Estimate, y = clean_var, fill = var_type)) +
  geom_bar(stat="identity") +
  scale_y_discrete(drop = FALSE) +

  theme_bw() +
  ylab("Variable") +
  xlab('Coefficient (log-odds)') +
  scale_fill_manual(name='Variable type', values=c('deepskyblue3', 'red3')) +
  geom_vline(xintercept = 0, linetype="dotted",
             color = "black") +
  geom_errorbar(aes(xmin = Estimate - Std..Error * 1.96,
                    xmax = Estimate + Std..Error * 1.96),
                width = 0.2, color = "black")

ggsave(plot = clusterbased_hyperg_coef_barplot,
       filename = "/mfs/gdouglas/scripts/ocean_mag_hgt/display_items/Main_clusterbased_hyperg_lm_coef.pdf",
       device = "pdf", width = 7, height = 4, units = "in", dpi=600)

