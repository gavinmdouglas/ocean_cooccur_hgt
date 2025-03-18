rm(list = ls(all.names = TRUE))

library(ggplot2)
library(reshape2)
library(ComplexHeatmap)
library(gridGraphics)
library(cowplot)

clusterbased_hits <- read.table('/mfs/gdouglas/projects/ocean_hgt_zenodo/putative_hgt/cluster/cross_level_tallies_norm.tsv.gz',
                                stringsAsFactors = FALSE, sep = '\t', header = TRUE)

combined <- data.frame(matrix(NA, nrow = 3, ncol = 6))
colnames(combined) <- clusterbased_hits$Level
rownames(combined) <- c("clusterbased_95", "clusterbased_99", "num_comparisons")

combined["clusterbased_95", ] <- clusterbased_hits$cluster_95_99
combined["clusterbased_99", ] <- clusterbased_hits$cluster_99
combined["num_comparisons", ] <- clusterbased_hits$num_comparisons

combined_percent <- combined
for (i in 1:nrow(combined_percent)) {
  combined_percent[i, ] <- (combined_percent[i, ] / sum(combined_percent[i, ])) * 100
}

combined_formatted <- format(combined, big.mark = ',', justify = 'none', trim = TRUE)

count_heatmap <- Heatmap(matrix = as.matrix(combined_percent),

                                   col = circlize::colorRamp2(c(0, 100), c('white', 'firebrick3')),

                                   heatmap_legend_param = list(title = "Percent"),

                                   show_heatmap_legend = TRUE,

                                   row_title_rot=0,
                                   row_labels = c(">= 95% and < 99%", ">= 99%", "Number of    \ngenome        \ncomparisons"),
                                   row_names_gp = gpar(col="grey30"),

                                   cluster_rows = FALSE,
                                   cluster_columns = FALSE,
                                   column_names_gp = gpar(col="grey30"),

                                   row_names_side = 'left',
                                   row_dend_side = 'right',
                                   column_names_rot = 0,
                                   column_names_centered = TRUE,

                                   cell_fun = function(j, i, x, y, width, height, fill) {
                                     if(! is.na(combined_formatted[i, j] > 0))
                                       grid.text(combined_formatted[i, j], x, y, gp = gpar(fontsize = 10), just = 'centre')
                                   })

count_heatmap <- grid.grabExpr(draw(column_title = "", count_heatmap))


# Panel b
cluster_norm_hits <- clusterbased_hits[, c("Level", grep(".norm", colnames(clusterbased_hits), value = TRUE))]
colnames(cluster_norm_hits)[1] <- "Inter.level"
cluster_norm_hits_long <- reshape2::melt(data = cluster_norm_hits, id.vars = "Inter.level")
cluster_norm_hits_long$variable <- as.character(cluster_norm_hits_long$variable)
cluster_norm_hits_long$variable[which(cluster_norm_hits_long$variable == "cluster_95_99.norm")] <- "Identity >= 95% and < 99%"
cluster_norm_hits_long$variable[which(cluster_norm_hits_long$variable == "cluster_99.norm")] <- "Identity >= 99%"
cluster_norm_hits_long$variable[which(cluster_norm_hits_long$variable == "both.norm")] <- "All hits (Identity >= 95%)"

norm_hits_combined <- cluster_norm_hits_long

tax_levels <- c("Genus", "Family","Order", "Class", "Phylum", "Domain")
norm_hits_combined$Inter.level <- factor(norm_hits_combined$Inter.level, levels = tax_levels)

norm_hits_combined$variable <- factor(norm_hits_combined$variable, levels = c("All hits (Identity >= 95%)", "Identity >= 95% and < 99%", "Identity >= 99%"))

colnames(norm_hits_combined)[which(colnames(norm_hits_combined) == "variable")] <- "identity_cutoff"
colnames(norm_hits_combined)[which(colnames(norm_hits_combined) == "value")] <- "num_hgt_by_num_compare"

mixed_color <- (col2rgb("orange") +  col2rgb("cornflowerblue")) / 2
mixed_color_hex <- rgb(mixed_color[1,]/255, mixed_color[2,]/255, mixed_color[3,]/255)

norm_hits_plot <- ggplot(data = norm_hits_combined,
                         aes(x = Inter.level, y = num_hgt_by_num_compare, colour = identity_cutoff)) +
  geom_point(size = 5, alpha=1) +
  scale_y_continuous(trans = "log10", labels = scales::label_number(), limits=c(0.0000005, 1)) +
  theme_bw() +
  scale_colour_manual(values = c(mixed_color_hex, "orange", "cornflowerblue")) +
  labs(colour = "Identity cut-off") +
  xlab("Inter-taxonomic level") +
  ylab(expression(atop(displaystyle(frac("No. putative HGT events", "No. genome comparisons")),
                       displaystyle("(log"[10] * " scale)")))) +
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5),
        legend.position = "inside",
        legend.position.inside = c(0.55, 0.85),
        legend.box = "horizontal",
        legend.background = element_rect(
          fill = "white",
          color = "gray80",
        ),
        legend.margin = margin(6, 6, 6, 6),
        legend.box.background = element_rect(
          fill = "white",
          color = NA,
        ),
        plot.title = element_text(size = 20, face = "bold"),
        plot.subtitle = element_text(size = 16),
        plot.caption = element_text(size = 10),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size=12)
  )


# Combine plot.
top_row <- plot_grid(NULL, count_heatmap, NULL, rel_widths = c(1, 5, 0.1), nrow=1, labels=c('', 'a', ''))
bottom_row <- plot_grid(norm_hits_plot, NULL, labels=c('b', ''), nrow=1, rel_widths = c(11, 1))
combined_plot <- plot_grid(top_row, bottom_row, nrow=2, rel_heights = c(1, 2))

ggsave(plot = combined_plot,
       filename = "/mfs/gdouglas/scripts/ocean_mag_hgt/display_items/Main_identity_tallies_and_prop.pdf",
       device = "pdf", width = 8.65, height = 8, units = "in", dpi=600)
