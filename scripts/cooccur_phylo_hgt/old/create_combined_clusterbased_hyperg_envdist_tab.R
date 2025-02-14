rm(list = ls(all.names = TRUE))

combined_file="/mfs/gdouglas/projects/ocean_mags/networks/allsamples/clusterbased/metaG_hyperg_cooccur.allsamples.combined.tsv.gz"
median_diff_file="/mfs/gdouglas/projects/ocean_mags/networks/combined_tables/allsamples_present_metadata_median_pairwise_diff.tsv.gz"

combined_info <- read.table(combined_file, header=TRUE, sep="\t", stringsAsFactors = FALSE, row.names = 1)
combined_info <- combined_info[which(! combined_info$diff_tax_level %in% c('Species', 'Strain')), ]
combined_info <- combined_info[which(rowSums(is.na(combined_info)) == 0), ]

combined_info$cooccur_ratio <- combined_info$cooccur_obs / combined_info$cooccur_exp
combined_info$cooccur <- 0
combined_info$cooccur[which(combined_info$cooccur_BH < 0.05 & combined_info$cooccur_ratio > 1)] <- 1

combined_info$hgt <- 0
combined_info$hgt[which(combined_info[, 'both_gene_count'] > 0)] <- 1

combined_info$tip_dist_orderedNorm <- bestNormalize::orderNorm(x = combined_info$tip_dist)$x.t

# Add in metadata as well.
pairwise_metadata_dist <- read.table(median_diff_file, header=TRUE, sep="\t", stringsAsFactors = FALSE, row.names = 1)

combined_info$depth <- pairwise_metadata_dist[rownames(combined_info), 'depth']
combined_info$latitude <- pairwise_metadata_dist[rownames(combined_info), 'latitude']
combined_info$longitude <- pairwise_metadata_dist[rownames(combined_info), 'longitude']
combined_info$temperature <- pairwise_metadata_dist[rownames(combined_info), 'temperature']
combined_info$oxygen <- pairwise_metadata_dist[rownames(combined_info), 'oxygen']
combined_info$salinity <- pairwise_metadata_dist[rownames(combined_info), 'salinity']

combined_info$depth_orderedNorm <-  bestNormalize::orderNorm(x = combined_info$depth)$x.t
combined_info$latitude_orderedNorm <-  bestNormalize::orderNorm(x = combined_info$latitude)$x.t
combined_info$longitude_orderedNorm <-  bestNormalize::orderNorm(x = combined_info$longitude)$x.t
combined_info$temperature_orderedNorm <-  bestNormalize::orderNorm(x = combined_info$temperature)$x.t
combined_info$oxygen_orderedNorm <-  bestNormalize::orderNorm(x = combined_info$oxygen)$x.t
combined_info$salinity_orderedNorm <-  bestNormalize::orderNorm(x = combined_info$salinity)$x.t

combined_info <- combined_info[which(rowSums(is.na(combined_info)) == 0), ]

write.table(x=combined_info,
            file='/mfs/gdouglas/projects/ocean_mags/networks/allsamples/clusterbased/metaG_hyperg_cooccur.allsamples.combined.filt.w_env_dist.tsv',
            quote=FALSE, sep = '\t', row.names=TRUE, col.names=NA)
