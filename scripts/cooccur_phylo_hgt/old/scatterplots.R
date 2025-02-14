rm(list = ls(all.names = TRUE))

library(ggplot2)
library(hexbin)
library(ggpointdensity)

simple_combined <- read.table("/mfs/gdouglas/projects/ocean_mags/coverm/network_working/metaG_simple_cooccur.combined.tsv.gz",
                              header = TRUE, sep = "\t", stringsAsFactors = FALSE)

simple_combined <- simple_combined[which(simple_combined$both_gene_count >= 1), ]
simple_combined <- simple_combined[which(! is.na(simple_combined$tip_dist)), ]
simple_combined <- simple_combined[which(! is.na(simple_combined$tip_dist)), ]
simple_combined[which(is.na(simple_combined$cooccur_simple_cooccur)), "cooccur_simple_cooccur"] <- 0

ggplot(data = simple_combined, aes(x = log(tip_dist + 0.0001), y = cooccur_simple_cooccur)) +
  stat_binhex() +
  geom_density_2d() +
  labs(x = "log(Tip dist. + 0.0001)", y = "Co-occurrence") +
  theme_bw()

ggplot(data = simple_combined, aes(x = log(tip_dist + 0.0001), y = cooccur_simple_cooccur)) +
  geom_pointdensity() +
  scale_colour_continuous(low="grey80", high="darkblue") +
  labs(x = "log(Tip dist. + 0.0001)", y = "Co-occurrence") +
  theme_bw()

simple_combined$both_gene_count_ceil <- simple_combined$both_gene_count
simple_combined$both_gene_count_ceil[which(simple_combined$both_gene_count_ceil >= 100)] <- 100
ggplot(data = simple_combined, aes(x = log(tip_dist + 0.0001), y = both_gene_count_ceil)) +
  geom_pointdensity() +
  scale_colour_continuous(low="grey80", high="darkblue") +
  labs(x = "log(Tip dist. + 0.0001)", y = "HGT [# genes]") +
  theme_bw()

ggplot(data = simple_combined, aes(x = cooccur_simple_cooccur, y = both_gene_count_ceil)) +
  geom_pointdensity() +
  scale_colour_continuous(low="grey80", high="darkblue") +
  labs(x = "Co-occurrence", y = "HGT (# genes)") +
  theme_bw()
