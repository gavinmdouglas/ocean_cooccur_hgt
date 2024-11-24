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
hyperg_combined_info$cooccur[which(hyperg_combined_info$cooccur_BH < 0.05 & hyperg_combined_info$cooccur_ratio > 1)] <- 'Yes'

hyperg_combined_info$hgt <- 'No'
hyperg_combined_info$hgt[which(hyperg_combined_info$both_gene_count > 0)] <- 'Yes'

hyperg_combined_info$hgt_relationship <- 'No HGT'
hyperg_combined_info$hgt_relationship[which(hyperg_combined_info$hgt == 'Yes')] <- 'HGT'
hyperg_combined_info$hgt_relationship <- factor(hyperg_combined_info$hgt_relationship, levels = c('No HGT', 'HGT'))

hyperg_combined_info <- hyperg_combined_info[which(rowSums(is.na(hyperg_combined_info)) == 0), ]

# To speed up plotting, randomly subset 100,000 points from 'No' HGT.
hyperg_combined_info_subsampled <- hyperg_combined_info[which(hyperg_combined_info$hgt == 'No'), ]
hyperg_combined_info_subsampled <- hyperg_combined_info_subsampled[sample(1:nrow(hyperg_combined_info_subsampled), 100000), ]
hyperg_combined_info_subsampled <- rbind(hyperg_combined_info_subsampled, hyperg_combined_info[which(hyperg_combined_info$hgt == 'Yes'), ])

pairwise_env_dist <- read.table('/mfs/gdouglas/projects/ocean_mags/networks/combined_tables/allsamples_present_metadata_median_pairwise_diff.tsv.gz',
                                header=TRUE, sep="\t", stringsAsFactors = FALSE, row.names = 1)

# Quick comparisons of orderedNorm vs. original.
# pairwise_env_dist$depth_orderedNorm <-  bestNormalize::orderNorm(x = pairwise_env_dist[, 'depth'])$x.t
# pairwise_env_dist$latitude_orderedNorm <-  bestNormalize::orderNorm(x = pairwise_env_dist[, 'latitude'])$x.t
# pairwise_env_dist$longitude_orderedNorm <-  bestNormalize::orderNorm(x = pairwise_env_dist[, 'longitude'])$x.t
# pairwise_env_dist$temperature_orderedNorm <-  bestNormalize::orderNorm(x = pairwise_env_dist[, 'temperature'])$x.t
# pairwise_env_dist$oxygen_orderedNorm <-  bestNormalize::orderNorm(x = pairwise_env_dist[, 'oxygen'])$x.t
# pairwise_env_dist$salinity_orderedNorm <-  bestNormalize::orderNorm(x = pairwise_env_dist[, 'salinity'])$x.t
# pairwise_env_dist_subsampled <- pairwise_env_dist[sample(1:nrow(pairwise_env_dist), 100000),]
#
par(mfrow=c(2,3))
plot(pairwise_env_dist_subsampled$depth, pairwise_env_dist_subsampled$depth_orderedNorm, main='Depth (pairwise difference)', ylab='orderedNorm', xlab='Original')
plot(pairwise_env_dist_subsampled$latitude, pairwise_env_dist_subsampled$latitude_orderedNorm, main='Latitude (pairwise difference)', ylab='orderedNorm', xlab='Original')
plot(pairwise_env_dist_subsampled$longitude, pairwise_env_dist_subsampled$longitude_orderedNorm, main='Longitude (pairwise difference)', ylab='orderedNorm', xlab='Original')
plot(pairwise_env_dist_subsampled$temperature, pairwise_env_dist_subsampled$temperature_orderedNorm, main='Temperature (pairwise difference)', ylab='orderedNorm', xlab='Original')
plot(pairwise_env_dist_subsampled$oxygen, pairwise_env_dist_subsampled$oxygen_orderedNorm, main='Oxygen (pairwise difference)', ylab='orderedNorm', xlab='Original')
plot(pairwise_env_dist_subsampled$salinity, pairwise_env_dist_subsampled$salinity_orderedNorm, main='Salinity (pairwise difference)', ylab='orderedNorm', xlab='Original')

prepped_tab <- data.frame(hgt_relationship=rep(hyperg_combined_info_subsampled$hgt_relationship, times = 6),
                          value=c(pairwise_env_dist[hyperg_combined_info_subsampled$taxa_combo, 'depth'],
                                  pairwise_env_dist[hyperg_combined_info_subsampled$taxa_combo, 'latitude'],
                                  pairwise_env_dist[hyperg_combined_info_subsampled$taxa_combo, 'longitude'],
                                  pairwise_env_dist[hyperg_combined_info_subsampled$taxa_combo, 'temperature'],
                                  pairwise_env_dist[hyperg_combined_info_subsampled$taxa_combo, 'oxygen'],
                                  pairwise_env_dist[hyperg_combined_info_subsampled$taxa_combo, 'salinity']),
                          variable=c(rep('Depth (diff.)', nrow(hyperg_combined_info_subsampled)),
                                     rep('Latitude (diff.)', nrow(hyperg_combined_info_subsampled)),
                                     rep('Longitude (diff.)', nrow(hyperg_combined_info_subsampled)),
                                     rep('Temperature (diff.)', nrow(hyperg_combined_info_subsampled)),
                                     rep('Oxygen (diff.)', nrow(hyperg_combined_info_subsampled)),
                                     rep('Salinity (diff.)', nrow(hyperg_combined_info_subsampled))))


env_by_cooccur_and_hgt <- ggplot(data = prepped_tab, aes(x = hgt_relationship, y = value, fill = hgt_relationship)) +
  geom_violin(fill='grey85', col='grey85') +
  geom_boxplot(outlier.shape=NA, alpha=0.3) +
  scale_fill_manual(name='HGT', values=c("grey80", "#56B4E9")) +
  facet_wrap(. ~ variable, scales='free_y') +
  theme_bw() +
  ylab('Median difference') +
  xlab('Horizontal gene transfer (HGT) relationship') +
  theme(legend.position = "none")


ggsave(plot = env_by_cooccur_and_hgt,
       filename = "/mfs/gdouglas/scripts/ocean_mag_hgt/display_items/Main_clusterbased_hgt_cooccur_envdist_overview.pdf",
       device = "pdf", width = 8, height = 5, units = "in", dpi=600)
