rm(list = ls(all.names = TRUE))

library(ggplot2)
library(cowplot)

results <- read.table("/mfs/gdouglas/projects/water_mags/blast_output/pairwise_tally_summary.tsv.gz",
                      header = TRUE, sep = "\t", stringsAsFactors = FALSE)

#results <- results[-which(results$Highest_tax_diff %in% c("Species", "Strain")), ]

results$both_gene_count[which(results$both_gene_count == 0)] <- NA
results$X95_gene_count[which(results$X95_gene_count == 0)] <- NA
results$X99_gene_count[which(results$X99_gene_count == 0)] <- NA
results$X95_hit_count[which(results$X95_hit_count == 0)] <- NA
results$X99_hit_count[which(results$X99_hit_count == 0)] <- NA

results$both_gene_count <- ceiling(results$both_gene_count)
results$X95_gene_count <- ceiling(results$X95_gene_count)
results$X99_gene_count <- ceiling(results$X99_gene_count)

results$Highest_tax_diff <- factor(results$Highest_tax_diff, levels = c("Strain", "Species", "Genus", "Family", "Order", "Class", "Phylum", "Domain"))

overall_num_hits <- ggplot(data = results, aes(x = Highest_tax_diff, y = both_hit_count)) +
  geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "orange") +
  scale_y_log10() +
  ylab(expression(log[10]*"(Number of BLAST hits)")) +
  xlab("Inter-taxonomic level") +
  theme_bw()

overall_gene_hits <- ggplot(data = results[which(! is.na(results$both_gene_count)), ],
                            aes(x = Highest_tax_diff, y = both_gene_count)) +
      geom_boxplot() +
      stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "cornflowerblue") +
      scale_y_log10() +
      ylab(expression(log[10]*"(Number of genes in hits)")) +
      xlab("Inter-taxonomic level") +
      theme_bw()

overall_plot_title <- ggdraw() + 
  draw_label("All BLAST hits (>= 95%)", x = 0, hjust = 0) +
  theme(plot.margin = margin(0, 0, 0, 320))

overall_plot <- plot_grid(overall_num_hits, overall_gene_hits)

overall_plot <- plot_grid(overall_plot_title, overall_plot,
                          rel_heights = c(0.1, 1), nrow = 2)


X95_num_hits <- ggplot(data = results[which(! is.na(results$X95_hit_count)), ],
                       aes(x = Highest_tax_diff, y = X95_hit_count)) +
  geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "orange") +
  scale_y_log10() +
  ylab(expression(log[10]*"(Number of BLAST hits)")) +
  xlab("Inter-taxonomic level") +
  theme_bw()

X95_gene_hits <- ggplot(data = results[which(! is.na(results$X95_gene_count)), ],
                        aes(x = Highest_tax_diff, y = X95_gene_count)) +
  geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "cornflowerblue") +
  scale_y_log10() +
  ylab(expression(log[10]*"(Number of genes in hits)")) +
  xlab("Inter-taxonomic level") +
  theme_bw()

X95_plot_title <- ggdraw() + 
  draw_label("Identity >= 95% and < 99%", x = 0, hjust = 0) +
  theme(plot.margin = margin(0, 0, 0, 320))

X95_plot <- plot_grid(X95_num_hits, X95_gene_hits)

X95_plot <- plot_grid(X95_plot_title, X95_plot,
                          rel_heights = c(0.1, 1), nrow = 2)


X99_num_hits <- ggplot(data = results[which(! is.na(results$X99_hit_count)), ],
                       aes(x = Highest_tax_diff, y = X99_hit_count)) +
  geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "orange") +
  scale_y_log10() +
  ylab(expression(log[10]*"(Number of BLAST hits)")) +
  xlab("Inter-taxonomic level") +
  theme_bw()

X99_gene_hits <- ggplot(data = results[which(! is.na(results$X99_gene_count)), ],
                        aes(x = Highest_tax_diff, y = X99_gene_count)) +
  geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "cornflowerblue") +
  scale_y_log10() +
  ylab(expression(log[10]*"(Number of genes in hits)")) +
  xlab("Inter-taxonomic level") +
  theme_bw()

X99_plot_title <- ggdraw() + 
  draw_label("Identity >= 99%", x = 0, hjust = 0) +
  theme(plot.margin = margin(0, 0, 0, 320))


X99_plot <- plot_grid(X99_num_hits, X99_gene_hits)

X99_plot <- plot_grid(X99_plot_title, X99_plot,
                      rel_heights = c(0.1, 1), nrow = 2)
