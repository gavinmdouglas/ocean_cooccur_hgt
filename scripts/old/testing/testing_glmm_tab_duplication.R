rm(list = ls(all.names = TRUE))

library(glmmTMB)
library(Matrix)
library(dplyr)

combined_file="/mfs/gdouglas/projects/ocean_mags/networks/allsamples/clusterbased/metaG_hyperg_cooccur.allsamples.combined.tsv.gz"
#combined_file="/mfs/gdouglas/projects/ocean_mags/networks/allsamples/clusterbased/metaG_hyperg_cooccur.allsamples.combined.filt.w_env_dist_head1000000.tsv"
cooccur_approach='hyperg'
hgt_tally_col='both_gene_count'
median_diff_file="/mfs/gdouglas/projects/ocean_mags/networks/combined_tables/allsamples_present_metadata_median_pairwise_diff.tsv.gz"
outfolder="/mfs/gdouglas/tmp/glmm_test/"

num_cores=8
keep_lower_levels = FALSE
suffix='.rds'

  if (! cooccur_approach %in% c('hyperg', 'simple', 'propr')) { stop('Co-occur approach must be hyperg, simple, or propr.') }
  if (length(grep(cooccur_approach, combined_file)) == 0) { stop('Co-occur approach not present in combined_file?') }

  combined_info <- read.table(combined_file, header=TRUE, sep="\t", stringsAsFactors = FALSE, row.names = 1)

  if (! keep_lower_levels) {
    combined_info <- combined_info[which(! combined_info$diff_tax_level %in% c('Species', 'Strain')), ]
  }
  combined_info <- combined_info[which(rowSums(is.na(combined_info)) == 0), ]

  if (cooccur_approach == 'hyperg') {
    combined_info$cooccur_ratio <- combined_info$cooccur_obs / combined_info$cooccur_exp
    combined_info$cooccur <- 0
    combined_info$cooccur[which(combined_info$cooccur_BH < 0.05 & combined_info$cooccur_ratio > 1)] <- 1
  } else if (cooccur_approach == 'simple') {
    combined_info <- combined_info[which(! is.na(combined_info$cooccur_simple_cooccur)), ]
    combined_info$cooccur <- bestNormalize::orderNorm(x = combined_info$cooccur_simple_cooccur)$x.t
  } else if (cooccur_approach == 'propr') {
    combined_info <- combined_info[which(! is.na(combined_info$cooccur_asso)), ]
    combined_info$cooccur <- bestNormalize::orderNorm(x = combined_info$cooccur_asso)$x.t
  }

  combined_info$hgt <- 0
  combined_info$hgt[which(combined_info[, hgt_tally_col] > 0)] <- 1
  glmm_family = "binomial"


  combined_info$tip_dist_orderedNorm <- bestNormalize::orderNorm(x = combined_info$tip_dist)$x.t

  # Add in metadata as well.
  pairwise_metadata_dist <- read.table(median_diff_file, header=TRUE, sep="\t", stringsAsFactors = FALSE, row.names = 1)

  combined_info$depth_orderedNorm <-  bestNormalize::orderNorm(x = pairwise_metadata_dist[rownames(combined_info), 'depth'])$x.t
  combined_info$latitude_orderedNorm <-  bestNormalize::orderNorm(x = pairwise_metadata_dist[rownames(combined_info), 'latitude'])$x.t
  combined_info$longitude_orderedNorm <-  bestNormalize::orderNorm(x = pairwise_metadata_dist[rownames(combined_info), 'longitude'])$x.t
  combined_info$temperature_orderedNorm <-  bestNormalize::orderNorm(x = pairwise_metadata_dist[rownames(combined_info), 'temperature'])$x.t
  combined_info$oxygen_orderedNorm <-  bestNormalize::orderNorm(x = pairwise_metadata_dist[rownames(combined_info), 'oxygen'])$x.t
  combined_info$salinity_orderedNorm <-  bestNormalize::orderNorm(x = pairwise_metadata_dist[rownames(combined_info), 'salinity'])$x.t

  combined_info <- combined_info[which(rowSums(is.na(combined_info)) == 0), ]

  hgt_by_cooccur <- glmmTMB(formula = hgt ~ cooccur,
                                        data = combined_info,
                                        family = glmm_family,
                                        control = glmmTMBControl(optimizer = nlminb,
                                                                 parallel = num_cores,
                                                                 optCtrl = list(iter.max = 300,
                                                                                eval.max = 400)))

  hgt_by_cooccur <- summary(hgt_by_cooccur)
  hgt_by_cooccur$call <- NULL
  saveRDS(object = hgt_by_cooccur, file = paste(outfolder, '/hgt_by_cooccur', suffix, sep = ''))
  rm(hgt_by_cooccur)

  col2keep <- c("hgt", "cooccur", "tip_dist_orderedNorm", "tip_dist_orderedNorm", "depth_orderedNorm", "latitude_orderedNorm",
                "longitude_orderedNorm", "temperature_orderedNorm", "oxygen_orderedNorm", "salinity_orderedNorm")
  tab_subset1 <- combined_info[, c(col2keep, 'taxon_i')]
  tab_subset2 <- combined_info[, c(col2keep, 'taxon_j')]

  colnames(tab_subset1) <- c(col2keep, 'taxon')
  colnames(tab_subset2) <- c(col2keep, 'taxon')

  tab_duplicated <- rbind(tab_subset1, tab_subset2)
  tab_duplicated$taxon <- factor(tab_duplicated$taxon)

  message('Ready to run first model')
  hgt_by_cooccur_w_dup_table_and_random <- glmmTMB(formula = hgt ~ cooccur + (1 | taxon),
                            data = tab_duplicated,
                            family = glmm_family,
                            control = glmmTMBControl(optimizer = nlminb,
                                                     parallel = num_cores,
                                                     optCtrl = list(iter.max = 300,
                                                                    eval.max = 400)))

  hgt_by_cooccur_w_dup_table_and_random <- summary(hgt_by_cooccur_w_dup_table_and_random)
  hgt_by_cooccur_w_dup_table_and_random$call <- NULL
  saveRDS(object = hgt_by_cooccur_w_dup_table_and_random, file = paste(outfolder, '/hgt_by_cooccur_w_dup_table_and_random', suffix, sep = ''))
  rm(hgt_by_cooccur_w_dup_table_and_random)

  hgt_by_cooccur_w_dup_table <- glmmTMB(formula = hgt ~ cooccur,
                                                   data = tab_duplicated,
                                                   family = glmm_family,
                                                   control = glmmTMBControl(optimizer = nlminb,
                                                                            parallel = num_cores,
                                                                            optCtrl = list(iter.max = 300,
                                                                                           eval.max = 400)))

  hgt_by_cooccur_w_dup_table <- summary(hgt_by_cooccur_w_dup_table)
  hgt_by_cooccur_w_dup_table$call <- NULL
  saveRDS(object = hgt_by_cooccur_w_dup_table, file = paste(outfolder, '/hgt_by_cooccur_w_dup_table', suffix, sep = ''))
  rm(hgt_by_cooccur_w_dup_table)
