library(glmmTMB)
library(lmerMultiMember)
library(Matrix)
library(dplyr)

# Function for computing GLMMs with growing number of fixed effects, which
# in the end was just a sanity check that nothing weird was going on when different
# numbers of fixed effects were added in.
# In addition, random effect of genome ID was dealt with by simply duplicating table and adding
# in genome ID for each paired observation once (as opposed to explicitly accounting for multi-membership,
# as below).
compute_hgt_cooccur_glmms <- function(combined_file,
                                      cooccur_approach,
                                      hgt_tally_col,
                                      median_diff_file,
                                      outfolder,
                                      num_cores=8,
                                      keep_lower_levels = FALSE,
                                      suffix='.rds') {

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

  if (cooccur_approach == 'simple') {
    combined_info$hgt <- bestNormalize::orderNorm(x = combined_info[, hgt_tally_col])$x.t
    glmm_family = "gaussian"
  } else {
    combined_info$hgt <- 0
    combined_info$hgt[which(combined_info[, hgt_tally_col] > 0)] <- 1
    glmm_family = "binomial"
  }

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

  col2keep <- c("hgt", "cooccur", "tip_dist_orderedNorm", "tip_dist_orderedNorm", "depth_orderedNorm", "latitude_orderedNorm",
                "longitude_orderedNorm", "temperature_orderedNorm", "oxygen_orderedNorm", "salinity_orderedNorm")
  tab_subset1 <- combined_info[, c(col2keep, 'taxon_i')]
  tab_subset2 <- combined_info[, c(col2keep, 'taxon_j')]

  colnames(tab_subset1) <- c(col2keep, 'taxon')
  colnames(tab_subset2) <- c(col2keep, 'taxon')
  tab_subset <- rbind(tab_subset1, tab_subset2)
  tab_subset$taxon <- factor(tab_subset$taxon)

  message('Ready to run first model')
  hgt_by_cooccur <- glmmTMB(formula = hgt ~ cooccur + (1 | taxon),
                            data = tab_subset,
                            family = glmm_family,
                            control = glmmTMBControl(optimizer = nlminb,
                                                     parallel = num_cores,
                                                     optCtrl = list(iter.max = 300,
                                                                    eval.max = 400)))

  hgt_by_cooccur <- summary(hgt_by_cooccur)
  hgt_by_cooccur$call <- NULL
  saveRDS(object = hgt_by_cooccur, file = paste(outfolder, '/hgt_by_cooccur', suffix, sep = ''))
  rm(hgt_by_cooccur)

  message('Ready to run second model')
  hgt_by_cooccur_and_tip <- glmmTMB(formula = hgt ~ cooccur + tip_dist_orderedNorm + (1 | taxon),
                            data = tab_subset,
                            family = glmm_family,
                            control = glmmTMBControl(optimizer = nlminb,
                                                     parallel = num_cores,
                                                     optCtrl = list(iter.max = 300,
                                                                    eval.max = 400)))
  hgt_by_cooccur_and_tip <- summary(hgt_by_cooccur_and_tip)
  hgt_by_cooccur_and_tip$call <- NULL
  saveRDS(object = hgt_by_cooccur_and_tip, file = paste(outfolder, '/hgt_by_cooccur_and_tip', suffix, sep = ''))
  rm(hgt_by_cooccur_and_tip)

  message('Ready to run third model')
  hgt_by_cooccur_tip_and_env <- glmmTMB(formula = hgt ~ cooccur + tip_dist_orderedNorm + depth_orderedNorm + latitude_orderedNorm + longitude_orderedNorm + temperature_orderedNorm + oxygen_orderedNorm + salinity_orderedNorm + (1 | taxon),
                                            data = tab_subset,
                                            family = glmm_family,
                                            control = glmmTMBControl(optimizer = nlminb,
                                                                     parallel = num_cores,
                                                                     optCtrl = list(iter.max = 300,
                                                                                    eval.max = 400)))
  hgt_by_cooccur_tip_and_env <- summary(hgt_by_cooccur_tip_and_env)
  hgt_by_cooccur_tip_and_env$call <- NULL
  saveRDS(object = hgt_by_cooccur_tip_and_env, file = paste(outfolder, '/hgt_by_cooccur_tip_and_env', suffix, sep = ''))

}



# Version to account for multi-membership genome IDs explicitly.
# Also, only fits model with all fixed effects.
compute_multimember_hgt_cooccur_glmms <- function(combined_file,
                                      cooccur_approach,
                                      hgt_tally_col,
                                      median_diff_file,
                                      outprefix,
                                      keep_lower_levels = FALSE) {

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

  if (cooccur_approach == 'simple') {
    combined_info$hgt <- bestNormalize::orderNorm(x = combined_info[, hgt_tally_col])$x.t
    glmm_family = "gaussian"
  } else {
    combined_info$hgt <- 0
    combined_info$hgt[which(combined_info[, hgt_tally_col] > 0)] <- 1
    glmm_family = "binomial"
  }

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

  col2keep <- c("hgt", "cooccur", "tip_dist_orderedNorm", "tip_dist_orderedNorm", "depth_orderedNorm", "latitude_orderedNorm",
                "longitude_orderedNorm", "temperature_orderedNorm", "oxygen_orderedNorm", "salinity_orderedNorm")

  # The following code is based on Ben Bolker's suggestion, discussed here: https://github.com/glmmTMB/glmmTMB/issues/1117
  # And this tutorial: https://jvparidon.github.io/lmerMultiMember/articles/lmermultimember_intro.html
  Wp <- lmerMultiMember::weights_from_columns(combined_info[, c("taxon_i", "taxon_j")])

  m1 <- lmerMultiMember::glmer(hgt ~ cooccur + tip_dist_orderedNorm + depth_orderedNorm + latitude_orderedNorm + longitude_orderedNorm + temperature_orderedNorm + oxygen_orderedNorm + salinity_orderedNorm + (1 | genome_pair),
                                family = "binomial",
                                memberships = list(genome_pair = Wp),
                                data = combined_info)

  m1_summary <- summary(m1)

  coefficients_out <- data.frame(m1_summary$coefficients)

  random_effects_out <- data.frame(broom.mixed::tidy(m1, effects = "ran_vals", conf.int = TRUE) %>%
                                          .[.$group == "genome_pair", ] %>%
                                          .[order(.$estimate, decreasing = TRUE), ])

  write.table(x=coefficients_out, file=paste(outprefix, '.coef.tsv', sep = ''), sep = '\t', row.names = TRUE, col.names = NA, quote=FALSE)
  write.table(x=random_effects_out, file=paste(outprefix, '.raneffects.tsv', sep = ''), sep = '\t', row.names = FALSE, col.names = TRUE, quote=FALSE)

}


compute_hgt_cooccur_w_subgroups <- function(combined_file,
                                            cooccur_approach,
                                            hgt_tally_col,
                                            median_diff_file,
                                            outfolder,
                                            num_cores=8,
                                            keep_lower_levels = FALSE,
                                            suffix='.rds',
                                            freeliving_genomes_file='/mfs/gdouglas/projects/ocean_hgt_zenodo/cooccur/genomes_freeliving_associated.tsv.gz',
                                            lessfiltered_genomes_file='/mfs/gdouglas/projects/ocean_hgt_zenodo/cooccur/genomes_lessfiltered_associated.tsv.gz',
                                            itermax=300,
                                            evalmax=400) {

  if (! cooccur_approach %in% c('hyperg', 'simple', 'propr')) { stop('Co-occur approach must be hyperg, simple, or propr.') }
  if (length(grep(cooccur_approach, combined_file)) == 0) { stop('Co-occur approach not present in combined_file?') }

  combined_info <- data.frame(data.table::fread(combined_file, sep = '\t', verbose = TRUE))
  rownames(combined_info) <- combined_info[, 1]
  combined_info <- combined_info[, -1]

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

  if (cooccur_approach == 'simple') {
    combined_info$hgt <- bestNormalize::orderNorm(x = combined_info[, hgt_tally_col])$x.t
    glmm_family = "gaussian"
  } else {
    combined_info$hgt <- 0
    combined_info$hgt[which(combined_info[, hgt_tally_col] > 0)] <- 1
    glmm_family = "binomial"
  }

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

  genomes_mainly_freeliving <- read.table(file = freeliving_genomes_file, header=FALSE, stringsAsFactors = FALSE)$V1
  genomes_mainly_lessfiltered <- read.table(file = lessfiltered_genomes_file, header=FALSE, stringsAsFactors = FALSE)$V1

  taxon_i_grouping <- rep('Other', nrow(combined_info))
  taxon_j_grouping <- rep('Other', nrow(combined_info))
  taxon_i_grouping[which(combined_info$taxon_i %in% genomes_mainly_freeliving)] <- 'freeliving'
  taxon_i_grouping[which(combined_info$taxon_i %in% genomes_mainly_lessfiltered)] <- 'lessfiltered'
  taxon_j_grouping[which(combined_info$taxon_j %in% genomes_mainly_freeliving)] <- 'freeliving'
  taxon_j_grouping[which(combined_info$taxon_j %in% genomes_mainly_lessfiltered)] <- 'lessfiltered'

  combined_info$filter_group_match <- 'Mixed'
  combined_info[which(taxon_i_grouping == 'freeliving' & taxon_j_grouping == 'freeliving'), 'filter_group_match'] <- 'Free-living'
  combined_info[which(taxon_i_grouping == 'lessfiltered' & taxon_j_grouping == 'lessfiltered'), 'filter_group_match'] <- 'Less-filtered'


  col2keep <- c("hgt", "cooccur", "tip_dist_orderedNorm", "filter_group_match", "depth_orderedNorm", "latitude_orderedNorm",
                "longitude_orderedNorm", "temperature_orderedNorm", "oxygen_orderedNorm", "salinity_orderedNorm")
  tab_subset1 <- combined_info[, c(col2keep, 'taxon_i')]
  tab_subset2 <- combined_info[, c(col2keep, 'taxon_j')]

  colnames(tab_subset1) <- c(col2keep, 'taxon')
  colnames(tab_subset2) <- c(col2keep, 'taxon')
  tab_subset <- rbind(tab_subset1, tab_subset2)
  tab_subset$taxon <- factor(tab_subset$taxon)

  hgt_by_cooccur_tip_and_env <- glmmTMB(formula = hgt ~ cooccur + tip_dist_orderedNorm + filter_group_match + depth_orderedNorm + latitude_orderedNorm + longitude_orderedNorm + temperature_orderedNorm + oxygen_orderedNorm + salinity_orderedNorm + (1 | taxon),
                                        data = tab_subset,
                                        family = glmm_family,
                                        control = glmmTMBControl(parallel = num_cores,
                                                                 optCtrl = list(iter.max = itermax,
                                                                                eval.max = evalmax)))

  hgt_assess_out <- performance::check_collinearity(hgt_by_cooccur_tip_and_env)
  print(hgt_assess_out)

  hgt_by_cooccur_tip_and_env <- summary(hgt_by_cooccur_tip_and_env)
  hgt_by_cooccur_tip_and_env$call <- NULL
  saveRDS(object = hgt_by_cooccur_tip_and_env, file = paste(outfolder, '/hgt_by_cooccur_tip_and_env_w_subgroups', suffix, sep = ''))

}
