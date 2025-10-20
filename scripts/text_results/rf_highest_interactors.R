rm(list = ls(all.names = TRUE))

pairwise_interactions <- read.table("/mfs/gdouglas/projects/ocean_hgt_zenodo/hgt_prev_analyses/random_forest_output/RF_sig_variable_pairwise_H_statistic.tsv.gz",
                                    sep = "\t", stringsAsFactors = FALSE, header=TRUE)

highest_var <- c('Oxygen', 'fCDOM', 'Nitrate')
subset_tab <- pairwise_interactions[pairwise_interactions$var1 %in% highest_var & pairwise_interactions$var2 %in% highest_var, ]
mean(subset_tab$interaction_strength)
sd(subset_tab$interaction_strength)