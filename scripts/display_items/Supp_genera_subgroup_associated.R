rm(list = ls(all.names = TRUE))

library(cowplot)
library(ggplot2)

lessfiltered_genomes <- read.table('~/projects/ocean_hgt_zenodo/cooccur/genomes_lessfiltered_associated.tsv.gz', stringsAsFactors = FALSE)$V1
freeliving_genomes <- read.table('~/projects/ocean_hgt_zenodo/cooccur/genomes_freeliving_associated.tsv.gz', stringsAsFactors = FALSE)$V1

taxonomy <- read.table("/mfs/gdouglas/projects/ocean_hgt_zenodo/mapfiles/MAG_taxa_breakdown.tsv.gz",
                       header = TRUE, stringsAsFactors = FALSE, sep = "\t", row.names=2)

lessfiltered_genera_tally <- table(taxonomy[lessfiltered_genomes, 'Genus'])
freeliving_genera_tally <- table(taxonomy[freeliving_genomes, 'Genus'])

lessfiltered_genera_tally <- data.frame(genus=names(lessfiltered_genera_tally), count=as.integer(lessfiltered_genera_tally))
lessfiltered_genera_tally <- lessfiltered_genera_tally[order(lessfiltered_genera_tally$count, decreasing=TRUE), ]
lessfiltered_genera_tally_top30 <- head(lessfiltered_genera_tally, 30)

freeliving_genera_tally <- data.frame(genus=names(freeliving_genera_tally), count=as.integer(freeliving_genera_tally))
freeliving_genera_tally <- freeliving_genera_tally[order(freeliving_genera_tally$count, decreasing=TRUE), ]
freeliving_genera_tally_top30 <- head(freeliving_genera_tally, 30)

lessfiltered_genera_tally_top30$genus <- factor(lessfiltered_genera_tally_top30$genus, levels = rev(lessfiltered_genera_tally_top30$genus))
lessfiltered_panel <- ggplot(data=lessfiltered_genera_tally_top30, aes(y = genus, x = count)) +
  geom_bar(stat="identity") +
  ylab('Genus') +
  xlab('Number of genomes') +
  theme_bw() +
  ggtitle('\"Less filtered\" sample associated') +
  theme(plot.title = element_text(hjust = 0.5))

freeliving_genera_tally_top30$genus <- factor(freeliving_genera_tally_top30$genus, levels = rev(freeliving_genera_tally_top30$genus))
freeliving_panel <- ggplot(data=freeliving_genera_tally_top30, aes(y = genus, x = count)) +
  geom_bar(stat="identity") +
  ylab('Genus') +
  xlab('Number of genomes') +
  theme_bw() +
  ggtitle('\"Free-living\" sample associated') +
  theme(plot.title = element_text(hjust = 0.5))

combined_plot <- plot_grid(freeliving_panel, lessfiltered_panel, labels=c('a', 'b'), nrow=2)

ggsave(filename='/mfs/gdouglas/scripts/ocean_mag_hgt/display_items/Supp_genera_subgroup_associated.pdf',
       plot = combined_plot,
       height = 9, width = 11, dpi=600, device="pdf")
