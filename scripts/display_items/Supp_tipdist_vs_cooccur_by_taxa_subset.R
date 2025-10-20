rm(list = ls(all.names = TRUE))

library(cowplot)
library(ggplot2)

hyperg_combined_info <- read.table("/mfs/gdouglas/projects/ocean_mags/networks/allsamples/clusterbased/metaG_hyperg_cooccur.allsamples.combined.tsv.gz",
                                   header=TRUE, sep="\t", stringsAsFactors = FALSE)
#hyperg_combined_info <- read.table("/mfs/gdouglas/projects/ocean_mags/networks/allsamples/clusterbased/metaG_hyperg_cooccur.allsamples.combined_head100000.tsv.gz",
#                                   header=TRUE, sep="\t", stringsAsFactors = FALSE)

hyperg_combined_info$cooccur_ratio <- hyperg_combined_info$cooccur_obs / hyperg_combined_info$cooccur_exp

hyperg_combined_info$cooccur_ratio_w_dummy <- (hyperg_combined_info$cooccur_obs + 1) / (hyperg_combined_info$cooccur_exp + 1)

hyperg_combined_info <- hyperg_combined_info[which(! is.na(hyperg_combined_info$tip_dist)), ]

hyperg_combined_info_higheronly <- hyperg_combined_info[which(! hyperg_combined_info$diff_tax_level %in% c('Species', 'Strain')), ]
spearman_higher_only <- cor.test(hyperg_combined_info_higheronly$cooccur_ratio, hyperg_combined_info_higheronly$tip_dist, method = "spearman")

hyperg_combined_info_loweronly <- hyperg_combined_info[which(hyperg_combined_info$diff_tax_level %in% c('Species', 'Strain')), ]
spearman_lower_only <- cor.test(hyperg_combined_info_loweronly$cooccur_ratio, hyperg_combined_info_loweronly$tip_dist, method = "spearman")


higheronly_dataset_hex <- ggplot(hyperg_combined_info_higheronly, aes(x=tip_dist, y=log2(cooccur_ratio_w_dummy))) +
  geom_hex(bins = 30) +
  scale_fill_viridis_c(name='Count') +
  theme_bw() +
  labs(title = paste0("Above genus only\n(Spearman's ρ = ", round(spearman_higher_only$estimate, 3), ", P < 0.001)"),
       y = expression("Co-occurrence ("*log[2] ~ "enrichment)"),
       x = 'Phylogenetic distance') +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylim(-7, 4.5)

loweronly_dataset_hex <- ggplot(hyperg_combined_info_loweronly, aes(x=tip_dist, y=log2(cooccur_ratio_w_dummy))) +
  geom_hex(bins = 30) +
  scale_fill_viridis_c(name='Count') +
  theme_bw() +
  labs(title = paste0("Within genus or species only\n(Spearman's ρ = ", round(spearman_lower_only$estimate, 3),", P < 0.001)"),
       y = expression("Co-occurrence ("*log[2] ~ "enrichment)"),
       x = 'Phylogenetic distance') +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylim(-7, 4.5)

combined_plot <- cowplot::plot_grid(higheronly_dataset_hex, loweronly_dataset_hex,
                                    labels = c('a', 'b'), nrow=1)

ggsave(plot = combined_plot,
       filename = "/mfs/gdouglas/scripts/ocean_mag_hgt/display_items/Supp_tipdist_vs_cooccur_by_taxa_subsets.png",
       device = "png", width = 10, height = 4, units = "in", dpi=400)
