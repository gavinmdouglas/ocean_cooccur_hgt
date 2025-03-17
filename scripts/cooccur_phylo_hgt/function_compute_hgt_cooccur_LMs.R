# Function for computing linear models, with different variables included.

library(glmmTMB)

compute_hgt_cooccur_LMs <- function(prepped_table_path,
                                    outfolder,
                                    num_cores=8) {

  intab <- read.table(prepped_table_path, header=TRUE, sep="\t", stringsAsFactors = FALSE)

  hgt_by_cooccur <- glm(hgt ~ cooccur,
                        data=intab, family=binomial(link="logit"))
  hgt_by_cooccur <- summary(hgt_by_cooccur)
  hgt_by_cooccur$call <- NULL
  hgt_by_cooccur$deviance.resid <- NULL
  saveRDS(object = hgt_by_cooccur, file = paste(outfolder, '/hgt_by_cooccur.glm.rds', sep = ''))

  hgt_by_cooccur_and_tip <- glm(hgt ~ cooccur + tip_dist_orderedNorm,
                                data=intab, family=binomial(link="logit"))
  hgt_by_cooccur_and_tip_assess <- performance::check_collinearity(hgt_by_cooccur_and_tip)
  write.table(x = data.frame(hgt_by_cooccur_and_tip_assess), file = paste(outfolder, '/hgt_by_cooccur_and_tip.glm.VIF.tsv', sep = ''), col.names = TRUE, quote = FALSE, sep = '\t', row.names = FALSE)
  hgt_by_cooccur_and_tip <- summary(hgt_by_cooccur_and_tip)
  hgt_by_cooccur_and_tip$call <- NULL
  hgt_by_cooccur_and_tip$deviance.resid <- NULL
  saveRDS(object = hgt_by_cooccur_and_tip, file = paste(outfolder, '/hgt_by_cooccur_and_tip.glm.rds', sep = ''))
  rm(hgt_by_cooccur_and_tip)

  hgt_by_cooccur_tip_and_env <- glm(hgt ~ cooccur + tip_dist_orderedNorm + depth_orderedNorm + latitude_orderedNorm + longitude_orderedNorm + temperature_orderedNorm + oxygen_orderedNorm + salinity_orderedNorm,
                                    data=intab, family=binomial(link="logit"))
  hgt_by_cooccur_tip_and_env_assess <- performance::check_collinearity(hgt_by_cooccur_tip_and_env)
  write.table(x = data.frame(hgt_by_cooccur_tip_and_env_assess), file = paste(outfolder, '/hgt_by_cooccur_tip_and_env.glm.VIF.tsv', sep = ''), col.names = TRUE, quote = FALSE, sep = '\t', row.names = FALSE)
  hgt_by_cooccur_tip_and_env <- summary(hgt_by_cooccur_tip_and_env)
  hgt_by_cooccur_tip_and_env$call <- NULL
  hgt_by_cooccur_tip_and_env$deviance.resid <- NULL
  saveRDS(object = hgt_by_cooccur_tip_and_env, file = paste(outfolder, '/hgt_by_cooccur_tip_and_env.glm.rds', sep = ''))

  # Also test effect of duplicating table, and add in genome number as random effect, as control, just to make sure the result is very similar.
  tab_subset1 <- intab
  tab_subset2 <- intab
  tab_subset1$taxon <- tab_subset1$taxon_i
  tab_subset2$taxon <- tab_subset1$taxon_j

  tab_dup <- rbind(tab_subset1, tab_subset2)
  tab_dup$taxon <- factor(tab_dup$taxon)
  glmm_hgt_by_cooccur_tip_and_env <- glmmTMB(formula = hgt ~ cooccur + tip_dist_orderedNorm + depth_orderedNorm + latitude_orderedNorm + longitude_orderedNorm + temperature_orderedNorm + oxygen_orderedNorm + salinity_orderedNorm,
                                        data = tab_dup,
                                        family = 'binomial',
                                        control = glmmTMBControl(parallel = num_cores, optCtrl = list(iter.max = 300, eval.max = 400)))
  glmm_hgt_by_cooccur_tip_and_env_assess <- performance::check_collinearity(glmm_hgt_by_cooccur_tip_and_env)
  write.table(x = data.frame(glmm_hgt_by_cooccur_tip_and_env_assess), file = paste(outfolder, '/hgt_by_cooccur_tip_and_env.glmm.VIF.tsv', sep = ''), col.names = TRUE, quote = FALSE, sep = '\t', row.names = FALSE)

  glmm_hgt_by_cooccur_tip_and_env <- summary(glmm_hgt_by_cooccur_tip_and_env)
  glmm_hgt_by_cooccur_tip_and_env$call <- NULL
  glmm_hgt_by_cooccur_tip_and_env$deviance.resid <- NULL
  saveRDS(object = glmm_hgt_by_cooccur_tip_and_env, file = paste(outfolder, '/hgt_by_cooccur_tip_and_env.glmm.rds', sep = ''))

  glmm_hgt_by_cooccur_tip_and_env_w_taxon <- glmmTMB(formula = hgt ~ cooccur + tip_dist_orderedNorm + depth_orderedNorm + latitude_orderedNorm + longitude_orderedNorm + temperature_orderedNorm + oxygen_orderedNorm + salinity_orderedNorm + (1 | taxon),
                                             data = tab_dup,
                                             family = 'binomial',
                                             control = glmmTMBControl(parallel = num_cores, optCtrl = list(iter.max = 300, eval.max = 400)))
  glmm_hgt_by_cooccur_tip_and_env_w_taxon_assess <- performance::check_collinearity(glmm_hgt_by_cooccur_tip_and_env_w_taxon)
  write.table(x = data.frame(glmm_hgt_by_cooccur_tip_and_env_w_taxon_assess), file = paste(outfolder, '/hgt_by_cooccur_tip_and_env_w_taxon.glmm.VIF.tsv', sep = ''), col.names = TRUE, quote = FALSE, sep = '\t', row.names = FALSE)
  glmm_hgt_by_cooccur_tip_and_env_w_taxon <- summary(glmm_hgt_by_cooccur_tip_and_env_w_taxon)
  glmm_hgt_by_cooccur_tip_and_env_w_taxon$call <- NULL
  glmm_hgt_by_cooccur_tip_and_env_w_taxon$deviance.resid <- NULL
  saveRDS(object = glmm_hgt_by_cooccur_tip_and_env_w_taxon, file = paste(outfolder, '/hgt_by_cooccur_tip_and_env_w_taxon.glmm.rds', sep = ''))

}

compute_hgt_cooccur_LMs_w_filtergroups <- function(prepped_table_path,
                                    outfolder,
                                    num_cores=8,
                                    freeliving_genomes_file,
                                    lessfiltered_genomes_file) {

  intab <- read.table(prepped_table_path, header=TRUE, sep="\t", stringsAsFactors = FALSE)

  genomes_mainly_freeliving <- read.table(file = freeliving_genomes_file, header=FALSE, stringsAsFactors = FALSE)$V1
  genomes_mainly_lessfiltered <- read.table(file = lessfiltered_genomes_file, header=FALSE, stringsAsFactors = FALSE)$V1

  taxon_i_grouping <- rep('Other', nrow(intab))
  taxon_j_grouping <- rep('Other', nrow(intab))
  taxon_i_grouping[which(intab$taxon_i %in% genomes_mainly_freeliving)] <- 'freeliving'
  taxon_i_grouping[which(intab$taxon_i %in% genomes_mainly_lessfiltered)] <- 'lessfiltered'
  taxon_j_grouping[which(intab$taxon_j %in% genomes_mainly_freeliving)] <- 'freeliving'
  taxon_j_grouping[which(intab$taxon_j %in% genomes_mainly_lessfiltered)] <- 'lessfiltered'

  intab$filter_group_match <- 'Mixed'
  intab[which(taxon_i_grouping == 'freeliving' & taxon_j_grouping == 'freeliving'), 'filter_group_match'] <- 'Free-living'
  intab[which(taxon_i_grouping == 'lessfiltered' & taxon_j_grouping == 'lessfiltered'), 'filter_group_match'] <- 'Less-filtered'

  intab$filter_group_match <- factor(intab$filter_group_match, levels = c('Mixed', 'Free-living', 'Less-filtered'))

  hgt_by_cooccur_and_tip <- glm(hgt ~ cooccur + tip_dist_orderedNorm + filter_group_match,
                                data=intab, family=binomial(link="logit"))
  hgt_by_cooccur_and_tip_assess <- performance::check_collinearity(hgt_by_cooccur_and_tip)
  write.table(x = data.frame(hgt_by_cooccur_and_tip_assess), file = paste(outfolder, '/hgt_by_cooccur_and_tip_w_filtergroup.glm.VIF.tsv', sep = ''), col.names = TRUE, quote = FALSE, sep = '\t', row.names = FALSE)
  hgt_by_cooccur_and_tip <- summary(hgt_by_cooccur_and_tip)
  hgt_by_cooccur_and_tip$call <- NULL
  hgt_by_cooccur_and_tip$deviance.resid <- NULL
  saveRDS(object = hgt_by_cooccur_and_tip, file = paste(outfolder, '/hgt_by_cooccur_and_tip_w_filtergroup.glm.rds', sep = ''))
  rm(hgt_by_cooccur_and_tip)

  hgt_by_cooccur_tip_and_env <- glm(hgt ~ cooccur + tip_dist_orderedNorm + filter_group_match + depth_orderedNorm + latitude_orderedNorm + longitude_orderedNorm + temperature_orderedNorm + oxygen_orderedNorm + salinity_orderedNorm,
                                    data=intab, family=binomial(link="logit"))
  hgt_by_cooccur_tip_and_env_assess <- performance::check_collinearity(hgt_by_cooccur_tip_and_env)
  write.table(x = data.frame(hgt_by_cooccur_tip_and_env_assess), file = paste(outfolder, '/hgt_by_cooccur_tip_and_env_w_filtergroup.glm.VIF.tsv', sep = ''), col.names = TRUE, quote = FALSE, sep = '\t', row.names = FALSE)
  hgt_by_cooccur_tip_and_env <- summary(hgt_by_cooccur_tip_and_env)
  hgt_by_cooccur_tip_and_env$call <- NULL
  hgt_by_cooccur_tip_and_env$deviance.resid <- NULL
  saveRDS(object = hgt_by_cooccur_tip_and_env, file = paste(outfolder, '/hgt_by_cooccur_tip_and_env_w_filtergroup.glm.rds', sep = ''))

  # Also test effect of duplicating table, and add in genome number as random effect, as control, just to make sure the result is very similar.
  tab_subset1 <- intab
  tab_subset2 <- intab
  tab_subset1$taxon <- tab_subset1$taxon_i
  tab_subset2$taxon <- tab_subset1$taxon_j

  tab_dup <- rbind(tab_subset1, tab_subset2)
  tab_dup$taxon <- factor(tab_dup$taxon)
  glmm_hgt_by_cooccur_tip_and_env <- glmmTMB(formula = hgt ~ cooccur + tip_dist_orderedNorm + filter_group_match + depth_orderedNorm + latitude_orderedNorm + longitude_orderedNorm + temperature_orderedNorm + oxygen_orderedNorm + salinity_orderedNorm,
                                             data = tab_dup,
                                             family = 'binomial',
                                             control = glmmTMBControl(parallel = num_cores, optCtrl = list(iter.max = 300, eval.max = 400)))
  glmm_hgt_by_cooccur_tip_and_env_assess <- performance::check_collinearity(glmm_hgt_by_cooccur_tip_and_env)
  write.table(x = data.frame(glmm_hgt_by_cooccur_tip_and_env_assess), file = paste(outfolder, '/hgt_by_cooccur_tip_and_env_w_filtergroup.glmm.VIF.tsv', sep = ''), col.names = TRUE, quote = FALSE, sep = '\t', row.names = FALSE)

  glmm_hgt_by_cooccur_tip_and_env <- summary(glmm_hgt_by_cooccur_tip_and_env)
  glmm_hgt_by_cooccur_tip_and_env$call <- NULL
  glmm_hgt_by_cooccur_tip_and_env$deviance.resid <- NULL
  saveRDS(object = glmm_hgt_by_cooccur_tip_and_env, file = paste(outfolder, '/hgt_by_cooccur_tip_and_env_w_filtergroup.glmm.rds', sep = ''))

  glmm_hgt_by_cooccur_tip_and_env_w_taxon <- glmmTMB(formula = hgt ~ cooccur + tip_dist_orderedNorm + filter_group_match + depth_orderedNorm + latitude_orderedNorm + longitude_orderedNorm + temperature_orderedNorm + oxygen_orderedNorm + salinity_orderedNorm + (1 | taxon),
                                                     data = tab_dup,
                                                     family = 'binomial',
                                                     control = glmmTMBControl(parallel = num_cores, optCtrl = list(iter.max = 300, eval.max = 400)))
  glmm_hgt_by_cooccur_tip_and_env_w_taxon_assess <- performance::check_collinearity(glmm_hgt_by_cooccur_tip_and_env_w_taxon)
  write.table(x = data.frame(glmm_hgt_by_cooccur_tip_and_env_w_taxon_assess), file = paste(outfolder, '/hgt_by_cooccur_tip_and_env_w_taxon_w_filtergroup.glmm.VIF.tsv', sep = ''), col.names = TRUE, quote = FALSE, sep = '\t', row.names = FALSE)
  glmm_hgt_by_cooccur_tip_and_env_w_taxon <- summary(glmm_hgt_by_cooccur_tip_and_env_w_taxon)
  glmm_hgt_by_cooccur_tip_and_env_w_taxon$call <- NULL
  glmm_hgt_by_cooccur_tip_and_env_w_taxon$deviance.resid <- NULL
  saveRDS(object = glmm_hgt_by_cooccur_tip_and_env_w_taxon, file = paste(outfolder, '/hgt_by_cooccur_tip_and_env_w_taxon_w_filtergroup.glmm.rds', sep = ''))

}
