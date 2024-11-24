rm(list = ls(all.names = TRUE))

library("glmmTMB")
library("performance")

intab <- data.frame(data.table::fread('/mfs/gdouglas/projects/ocean_mags/networks/allsamples/clusterbased/metaG_hyperg_cooccur.allsamples.combined.filt.w_env_dist.tsv',
                                     sep = '\t', verbose = TRUE))

# Run basic generalized linear model (without random effects),
# simply to get estimates of variance.

# First, get logistic model with response of co-occurrence vs. all env dist variables.
# Get variance inflation factors estimated, to help determine whether variables are sufficiently independent (or whether a decomposition approach for instance should be used).
# logistic_cooccur_by_envdist <- glmmTMB(formula = cooccur ~ depth_orderedNorm + latitude_orderedNorm + longitude_orderedNorm + temperature_orderedNorm + oxygen_orderedNorm + salinity_orderedNorm,
#                                        data = intab,
#                                        family = binomial(link = "logit"))
#
# assess_out <- performance::check_collinearity(logistic_cooccur_by_envdist)
# assess_out
# summary(assess_out)

genomes_mainly_freeliving <- read.table(file = '/mfs/gdouglas/projects/ocean_hgt_zenodo/cooccur/genomes_freeliving_associated.tsv.gz', header=FALSE, stringsAsFactors = FALSE)$V1
genomes_mainly_lessfiltered <- read.table(file = '/mfs/gdouglas/projects/ocean_hgt_zenodo/cooccur/genomes_lessfiltered_associated.tsv.gz', header=FALSE, stringsAsFactors = FALSE)$V1

taxon_i_grouping <- rep('Other', nrow(intab))
taxon_j_grouping <- rep('Other', nrow(intab))
taxon_i_grouping[which(intab$taxon_i %in% genomes_mainly_freeliving)] <- 'freeliving'
taxon_i_grouping[which(intab$taxon_i %in% genomes_mainly_lessfiltered)] <- 'lessfiltered'
taxon_j_grouping[which(intab$taxon_j %in% genomes_mainly_freeliving)] <- 'freeliving'
taxon_j_grouping[which(intab$taxon_j %in% genomes_mainly_lessfiltered)] <- 'lessfiltered'

intab$filter_group_match <- 'Mixed'
intab[which(taxon_i_grouping == 'freeliving' & taxon_j_grouping == 'freeliving'), 'filter_group_match'] <- 'Free-living'
intab[which(taxon_i_grouping == 'lessfiltered' & taxon_j_grouping == 'lessfiltered'), 'filter_group_match'] <- 'Less-filtered'

col2keep <- c("hgt", "cooccur", "tip_dist_orderedNorm", "tip_dist_orderedNorm", "depth_orderedNorm", "latitude_orderedNorm",
              "longitude_orderedNorm", "temperature_orderedNorm", "oxygen_orderedNorm", "salinity_orderedNorm", "filter_group_match")
tab_subset1 <- intab[, c(col2keep, 'taxon_i')]
tab_subset2 <- intab[, c(col2keep, 'taxon_j')]

colnames(tab_subset1) <- c(col2keep, 'taxon')
colnames(tab_subset2) <- c(col2keep, 'taxon')
tab_dup <- rbind(tab_subset1, tab_subset2)
tab_dup$taxon <- factor(tab_dup$taxon)

tab_dup$filter_group_match <- factor(tab_dup$filter_group_match, levels=c('Mixed', 'Free-living', 'Less-filtered'))

glmm_hgt_model <- glmmTMB(formula = hgt ~ cooccur + tip_dist_orderedNorm + filter_group_match + depth_orderedNorm + latitude_orderedNorm + longitude_orderedNorm + temperature_orderedNorm + oxygen_orderedNorm + salinity_orderedNorm + (1 | taxon),
                          data = tab_dup,
                          family = binomial(link = "logit"))

hgt_assess_out <- performance::check_collinearity(glmm_hgt_model)
hgt_assess_out

glmm_hgt_model_summary <- summary(glmm_hgt_model)
glmm_hgt_model_summary$call <- NULL

saveRDS(object = glmm_hgt_model_summary, file = '/mfs/gdouglas/projects/ocean_mags/glmm_working/clusterbased/hyperg/hgt_by_cooccur_tip_and_env_w_samplegroup.rds')
saveRDS(object = glmm_hgt_model, file = '/mfs/gdouglas/projects/ocean_mags/glmm_working/clusterbased/hyperg/hgt_by_cooccur_tip_and_env_w_samplegroup_FULL.rds')



# RUN multi-member test, rather than simplistic duplication of table.
# The following code is based on Ben Bolker's suggestion, discussed here: https://github.com/glmmTMB/glmmTMB/issues/1117
# And this tutorial: https://jvparidon.github.io/lmerMultiMember/articles/lmermultimember_intro.html

# First, try running on random subset of table, as running on full table may be prohibitively slow...

set.seed(102)
intab_subsample <- intab[sample(1:nrow(intab), 1000000), ]

Wp_subsample <- lmerMultiMember::weights_from_columns(intab_subsample[, c("taxon_i", "taxon_j")])

m1_subsample <- lmerMultiMember::glmer(hgt ~ cooccur + tip_dist_orderedNorm + depth_orderedNorm + filter_group_match + latitude_orderedNorm + longitude_orderedNorm + temperature_orderedNorm + oxygen_orderedNorm + salinity_orderedNorm + (1 | genome_pair),
                             family = "binomial",
                             memberships = list(genome_pair = Wp_subsample),
                             data = intab_subsample)

m1_subsample_summary <- summary(m1_subsample)

coefficients_out <- data.frame(m1_subsample_summary$coefficients)

random_effects_out <- data.frame(broom.mixed::tidy(m1_subsample, effects = "ran_vals", conf.int = TRUE) %>%
                                   .[.$group == "genome_pair", ] %>%
                                   .[order(.$estimate, decreasing = TRUE), ])

write.table(x=coefficients_out,
            file="/mfs/gdouglas/projects/ocean_mags/glmm_working/multimember/out_clusterbased_hyperg_w_subgroups_subsample1mill.coef.tsv",
            sep = '\t', row.names = TRUE, col.names = NA, quote=FALSE)

write.table(x=random_effects_out,
            file="/mfs/gdouglas/projects/ocean_mags/glmm_working/multimember/out_clusterbased_hyperg_w_subgroups_subsample1mill.raneffects.tsv",
            sep = '\t', row.names = TRUE, col.names = NA, quote=FALSE)


write.table(x=, file=paste(outprefix, '.raneffects.tsv', sep = ''), sep = '\t', row.names = FALSE, col.names = TRUE, quote=FALSE)


# col_of_interest <- c('cooccur_ratio', 'depth', 'latitude', 'longitude', 'temperature', 'oxygen', 'salinity')
#
# env_median_dist_vars <- c('depth', 'latitude', 'longitude', 'temperature', 'oxygen', 'salinity')
#
# residuals_list <- list()
# for(var in env_vars) {
#   model <- glm(as.formula(paste0("cooccur_ratio ~", var)),
#                family = binomial(link = "logit"),
#                data = intab)
#   residuals_list[[var]] <- residuals(model, type = "deviance")
# }
#
#
# cor_matrix <- cor(head(intab[, col_of_interest], 100),
#                   method = "spearman",
#                   use = "pairwise.complete.obs")

