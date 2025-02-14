rm(list = ls(all.names = TRUE))

model_info1 <- readRDS(file = '/mfs/gdouglas/projects/ocean_mags/glmm_working/tara_combined_summary.rds')
model_info2 <- readRDS(file = '/mfs/gdouglas/projects/ocean_mags/glmm_working/progenomes_combined_summary.rds')
model_info <- c(model_info1, model_info2)

# Make quick AIC table.
AIC_tab <- data.frame(matrix(NA, nrow = length(model_info), ncol = 5))
colnames(AIC_tab) <- c('genome_set', 'tool_and_subset', 'approach', 'rds', 'aic')

for (i in 1:length(model_info)) {
  AIC_tab[i, c(1, 2, 3, 4)] <- c(model_info[[i]]$genomes, model_info[[i]]$tool_and_subset, model_info[[i]]$approach, basename(names(model_info)[i]))
  AIC_tab[i, 5] <- model_info[[i]]$aic
}

AIC_tab[grep('clusterbased_geotraces_samples', AIC_tab$tool_and_subset), 'approach'] <- 'hyperg'
AIC_tab[grep('clusterbased_tara_samples', AIC_tab$tool_and_subset), 'approach'] <- 'hyperg'
AIC_tab[grep('clusterbased_hyperg', AIC_tab$tool_and_subset), 'approach'] <- 'hyperg'

AIC_tab_lm <- AIC_tab[grep('.glmm.rds$', AIC_tab$rds, invert = TRUE), ]

unique_data_sets <- AIC_tab_lm[, 1:3]
unique_data_sets <- unique_data_sets[order(unique_data_sets$genome_set, unique_data_sets$tool_and_subset, unique_data_sets$approach), ]
unique_data_sets <- unique_data_sets[-which(duplicated(unique_data_sets)), ]

best_model <- character()
delta_aic <- numeric()

for (j in 1:nrow(unique_data_sets)) {

  genome_set_j <- unique_data_sets[j, 'genome_set']
  tool_and_subset_j <- unique_data_sets[j, 'tool_and_subset']
  approach_j <- unique_data_sets[j, 'approach']

  combo_subset <- AIC_tab_lm[which(AIC_tab_lm$genome_set == genome_set_j & AIC_tab_lm$tool_and_subset == tool_and_subset_j & AIC_tab_lm$approach == approach_j), ]

  combo_subset <- combo_subset[order(combo_subset$aic, decreasing = FALSE), ]

  best_model <- c(best_model, gsub('.glm.rds', '', combo_subset$rds[1]))

  delta_aic <- c(delta_aic, combo_subset$aic[1] - combo_subset$aic[2])

}

unique_data_sets$best_model <- best_model
unique_data_sets$delta_aic <- delta_aic

unique_data_sets
