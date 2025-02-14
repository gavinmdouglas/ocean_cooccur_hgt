rm(list = ls(all.names = TRUE))

library(ggplot2)
library(cowplot)

read_in_glmm_summary_rds_and_plot <- function(filepath, title, num_digits=2) {
  variable_map <- list()
  variable_map[["(Intercept)"]] <- 'Intercept'
  variable_map[["cooccur"]] <- 'Co-occurrence'
  variable_map[["tip_dist_orderedNorm"]] <- 'Inter-tip distance'
  variable_map[["filter_group_matchFree-living"]] <- '\"Free-living\" enriched'
  variable_map[["filter_group_matchLess-filtered"]] <- '\"Less filtered\" enriched'
  variable_map[["filter_group_matchMixed"]] <- '\"Mixed\" enriched'
  variable_map[["depth_orderedNorm"]] <- 'Depth (diff.)'
  variable_map[["latitude_orderedNorm"]] <- 'Latitude (diff.)'
  variable_map[["longitude_orderedNorm"]] <- 'Longitude (diff.)'
  variable_map[["temperature_orderedNorm"]] <- 'Temperature (diff.)'
  variable_map[["oxygen_orderedNorm"]] <- 'Oxygen (diff.)'
  variable_map[["salinity_orderedNorm"]] <- 'Salinity (diff.)'

  var_order <- c('Co-occurrence', '\"Free-living\" enriched', '\"Less filtered\" enriched',
                 'Inter-tip distance', 'Depth (diff.)', 'Latitude (diff.)',
                 'Longitude (diff.)', 'Temperature (diff.)', 'Oxygen (diff.)', 'Salinity (diff.)')

  coef_tab <- data.frame(readRDS(filepath)$coefficients$cond)
  coef_tab$BH <- p.adjust(coef_tab$Pr...z.., 'BH')
  coef_tab$variable <- rownames(coef_tab)
  coef_tab$clean_var <- NA
  for (i in 1:nrow(coef_tab)) {
    coef_tab[i, 'clean_var'] <- variable_map[[coef_tab$variable[i]]]
  }

  intercept_tmp <- format(round(coef_tab[which(coef_tab$clean_var == 'Intercept'), 'Estimate'], num_digits), nsmall = 2)
  coef_tab <- coef_tab[which(coef_tab$clean_var != 'Intercept'), ]

  coef_tab$clean_var <- factor(coef_tab$clean_var, levels=rev(coef_tab$clean_var))

  plot_out <- ggplot(data = coef_tab, aes(x = Estimate, y = clean_var)) +
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
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5))

}

plot_simple_cluster_coef <- read_in_glmm_summary_rds_and_plot('/mfs/gdouglas/projects/ocean_mags/glmm_working/clusterbased/simple/hgt_by_cooccur_tip_and_env_w_subgroups.rds',
                                                              'Simple - cluster-based', 4)
plot_propr_cluster_coef <- read_in_glmm_summary_rds_and_plot('/mfs/gdouglas/projects/ocean_mags/glmm_working/clusterbased/propr/hgt_by_cooccur_tip_and_env_w_subgroups.rds',
                                                             'propr - cluster-based')
plot_hyperg_blast_coef <- read_in_glmm_summary_rds_and_plot('/mfs/gdouglas/projects/ocean_mags/glmm_working/blast/hyperg/hgt_by_cooccur_tip_and_env_w_subgroups.rds',
                                                            'HyperG - BLASTn')
plot_hyperg_ranger_coef <- read_in_glmm_summary_rds_and_plot('/mfs/gdouglas/projects/ocean_mags/glmm_working/rangerdtl/hyperg/hgt_by_cooccur_tip_and_env_w_subgroups.rds',
                                                             'HyperG - RANGER-DTL')

combined_plot <- plot_grid(plot_simple_cluster_coef, plot_propr_cluster_coef, plot_hyperg_blast_coef, plot_hyperg_ranger_coef,
                           labels=c('a', 'b', 'c', 'd'))

ggsave(plot = combined_plot,
       filename = "/mfs/gdouglas/scripts/ocean_mag_hgt/display_items/Supp_alternative_workflow_glmms.pdf",
       device = "pdf", width = 8, height = 6, units = "in", dpi=600)

