rm(list = ls(all.names = TRUE))

library(ggplot2)
library(reshape2)
library(ComplexHeatmap)
library(gridGraphics)
library(cowplot)

hyperg_combined_info <- read.table("/mfs/gdouglas/projects/ocean_mags/networks/allsamples/clusterbased/metaG_hyperg_cooccur.allsamples.combined.tsv.gz",
                                   header=TRUE, sep="\t", stringsAsFactors = FALSE)
hyperg_combined_info <- hyperg_combined_info[which(! hyperg_combined_info$diff_tax_level %in% c('Species', 'Strain')), ]

hyperg_combined_info$cooccur <- 'No'
hyperg_combined_info$cooccur_ratio <- hyperg_combined_info$cooccur_obs / hyperg_combined_info$cooccur_exp
hyperg_combined_info$cooccur[which(hyperg_combined_info$cooccur_BH < 0.05 & hyperg_combined_info$cooccur_ratio > 1)] <- 'Yes'

hyperg_combined_info$hgt <- 'No'
hyperg_combined_info$hgt[which(hyperg_combined_info$both_gene_count > 0)] <- 'Yes'

contingency_tab <- data.frame(matrix(NA, nrow=2, ncol=2))
colnames(contingency_tab) <- c('hgt_no', 'hgt_yes')
rownames(contingency_tab) <- c('cooccur_no', 'cooccur_yes')

contingency_tab['cooccur_no', 'hgt_no'] <- length(which(hyperg_combined_info$hgt == 'No' & hyperg_combined_info$cooccur == 'No'))
contingency_tab['cooccur_yes', 'hgt_no'] <- length(which(hyperg_combined_info$hgt == 'No' & hyperg_combined_info$cooccur == 'Yes'))
contingency_tab['cooccur_no', 'hgt_yes'] <- length(which(hyperg_combined_info$hgt == 'Yes' & hyperg_combined_info$cooccur == 'No'))
contingency_tab['cooccur_yes', 'hgt_yes'] <- length(which(hyperg_combined_info$hgt == 'Yes' & hyperg_combined_info$cooccur == 'Yes'))

contingency_formatted <- format(contingency_tab, big.mark = ',', justify = 'none', trim = TRUE)

contingency_percent <- contingency_tab
contingency_percent$hgt_no <- contingency_tab$hgt_no / sum(contingency_tab$hgt_no)
contingency_percent$hgt_yes <- contingency_tab$hgt_yes / sum(contingency_tab$hgt_yes)
contingency_percent <- contingency_percent * 100

contingency_heatmap <- Heatmap(matrix = as.matrix(contingency_percent),

                               col = circlize::colorRamp2(c(0, 100), c('white', 'firebrick3')),

                               heatmap_legend_param = list(title = "Percent"),

                               show_heatmap_legend = TRUE,

                               row_gap = unit(5, "mm"),
                               row_title_rot=0,
                               row_title = 'Co-occurring',
                               row_labels = c('No', 'Yes'),

                               cluster_rows = FALSE,
                               cluster_columns = FALSE,
                               column_title = 'Horizontal gene transfer\nrelationship (Cluster-based)',
                               column_labels = c('No', 'Yes'),
                               column_split = c('No', 'Yes'),
                               column_gap = unit(10, "mm"),

                               row_names_side = 'left',
                               row_dend_side = 'right',
                               column_names_rot = 0,
                               column_names_centered = TRUE,

                               cell_fun = function(j, i, x, y, width, height, fill) {
                                 if(! is.na(contingency_formatted[i, j] > 0))
                                   grid.text(contingency_formatted[i, j], x, y, gp = gpar(fontsize = 10), just = 'centre')
                               })

contingency_heatmap <- grid.grabExpr(draw(column_title = "", contingency_heatmap))

hyperg_combined_info$hgt_relationship <- 'No HGT'
hyperg_combined_info$hgt_relationship[which(hyperg_combined_info$hgt == 'Yes')] <- 'HGT'
hyperg_combined_info$hgt_relationship <- factor(hyperg_combined_info$hgt_relationship, levels = c('No HGT', 'HGT'))

tip_dist_by_cooccur_and_hgt <- ggplot(data = subset_tab, aes(x = cooccur, y = tip_dist)) +
                                      geom_violin(fill='grey85', col='grey85') +
                                      geom_boxplot(outlier.shape=NA, alpha=0.3) +
                                      facet_wrap(. ~ hgt_relationship) +
                                      theme_bw() +
                                      ylab('Inter-tip distance') +
                                      xlab('Co-occurring')

top_row <- plot_grid(NULL, contingency_heatmap, NULL, labels=c('', 'a', ''), nrow=1, rel_widths = c(1, 2, 1))

combined_plot <- plot_grid(top_row, tip_dist_by_cooccur_and_hgt,
                          labels = c('', 'b'), nrow = 2, rel_heights = c(1, 2))

ggsave(plot = combined_plot,
       filename = "/mfs/gdouglas/scripts/ocean_mag_hgt/display_items/Main_clusterbased_hgt_cooccur_tipdist_overview.pdf",
       device = "pdf", width = 9, height = 7, units = "in", dpi=600)
