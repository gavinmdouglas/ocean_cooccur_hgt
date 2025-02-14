rm(list = ls(all.names = TRUE))

source('/mfs/gdouglas/scripts/ocean_mag_hgt/scripts/functions.R')

library(glmmTMB)

combined <- read.table('/mfs/gdouglas/projects/ocean_mags/summary_files/gene_info_and_hgt.tsv.gz',
                       header = TRUE, sep = '\t', stringsAsFactors = FALSE)
combined <- combined[which(rowSums(is.na(combined)) == 0), ]

# First restrict to Panaroo gene families present at least three times.
species_gf_tallies <- table(combined$panaroo_sp_gf)
combined <- combined[which(combined$panaroo_sp_gf %in% names(species_gf_tallies)[which(species_gf_tallies >= 3)]), ]

species_gf_tallies <- table(combined$panaroo_sp_gf)

combined_panaroosp_and_COG <- combined[, c('COG_category', 'ranger_species', 'panaroo_sp_gf')]
combined_panaroosp_and_COG <- combined_panaroosp_and_COG[which(! duplicated(combined_panaroosp_and_COG)), ]

species_gf_tallies_no_hgt <- table(combined[which(combined$ranger_hgt == 'No'), 'panaroo_sp_gf'])
species_gf_tallies_w_hgt <- table(combined[which(combined$ranger_hgt == 'Yes'), 'panaroo_sp_gf'])

unique_panaroo_sp_gf <- unique(combined$panaroo_sp_gf)
prop_hgt <- data.frame(matrix(NA, nrow=length(unique_panaroo_sp_gf), ncol = 4))
colnames(prop_hgt) <- c('gene_family', 'num_hgt', 'num_non_hgt', 'prop_hgt')

prop_hgt$num_hgt <- 0
prop_hgt$num_non_hgt <- 0

prop_hgt[names(species_gf_tallies_no_hgt), 'num_non_hgt'] <- as.integer(species_gf_tallies_no_hgt)
prop_hgt[names(species_gf_tallies_w_hgt), 'num_hgt'] <- as.integer(species_gf_tallies_w_hgt)

combined_panaroosp_and_COG <- combined[, c('COG_category', 'ranger_species', 'panaroo_sp_gf')]
combined_panaroosp_and_COG <- combined_panaroosp_and_COG[-which(duplicated(combined_panaroosp_and_COG)), ]
rownames(combined_panaroosp_and_COG) <- combined_panaroosp_and_COG$panaroo_sp_gf

prop_hgt$ranger_sp <- combined_panaroosp_and_COG[rownames(prop_hgt), ]

rownames(prop_hgt) <- unique_panaroo_sp_gf
for (gene_family in unique_panaroo_sp_gf) {
  print(gene_family)
  tab_subset <- combined[which(combined$panaroo_sp_gf == gene_family), ]
  prop_hgt[gene_family, 'prop_hgt'] <- length(which(tab_subset$panaroo_sp_gf == gene_family)) / nrow(tab_subset)
}

# Then convert gene family count to
combined_COG$panaroo_sp_gf_tallies <- as.numeric(panaroo_sp_tallies[combined_COG$panaroo_sp_gf])
order_norm_out <- bestNormalize::orderNorm(as.numeric(panaroo_sp_tallies[combined_COG$panaroo_sp_gf]))
combined_COG$orderedNorm_panaroo_sp_gf_tallies <- order_norm_out$x.t



# First run test on COG categories.
combined_COG <- combined
combined_COG <- combined_COG[which(combined_COG$COG_category != '-'), ]
combined_COG <- split_multi_category_rows(in_df = combined_COG, category_col = 'COG_category', num_cores = 20)
combined_COG <- combined_COG[which(combined_COG$COG_category != '-'), ]

expected_hgt <- read.table('/mfs/gdouglas/projects/ocean_hgt_zenodo/mapfiles/Dmitrijeva2024_COG_category_HGT_summary.tsv.gz',
                           header = TRUE, sep = '\t', stringsAsFactors = FALSE)
expected_enriched <- expected_hgt$COG_category[which(expected_hgt$Custom_summary == 'Enriched')]
expected_enriched <- expected_enriched[which(expected_enriched != "K (reg.)")]
expected_depleted <- expected_hgt$COG_category[which(expected_hgt$Custom_summary == 'Depleted')]
expected_mixed <- expected_hgt$COG_category[which(expected_hgt$Custom_summary == 'Mixed')]
expected_mixed[expected_mixed == "K (all)"] <- "K"

to_ignore <- c("A", "B", "Y", "Z")

combined_COG <- combined_COG[which(! combined_COG$COG_category %in% to_ignore), ]

combined_COG[which(combined_COG$COG_category %in% expected_mixed), "COG_category"] <- 'Mixed_expectation'
unique_COG_categories <- unique(combined_COG$COG_category)
COG_category_levels <- c("Mixed_expectation", unique_COG_categories[which(unique_COG_categories != "Mixed_expectation")])
combined_COG$COG_category <- factor(combined_COG$COG_category, levels = COG_category_levels)

RANGER_subset_col_i <- which(colnames(combined_COG) == "ranger_hgt")

RANGER_COG_glmm_out <- run_glmms_per_grouping(in_tab = combined_COG,
                               subset_col_i = RANGER_subset_col_i,
                               dependent_var_string = "COG_category + (1 | orderedNorm_panaroo_sp_gf_tallies) + (1 | ranger_species) + (1 | panaroo_sp_gf)",
                               backup1_dependent_var_string = "COG_category + (1 | orderedNorm_panaroo_sp_gf_tallies) + (1 | ranger_species)",
                               backup2_dependent_var_string = "COG_category + (1 | orderedNorm_panaroo_sp_gf_tallies)",
                               convert_to_integer_no_yes=TRUE,
                               num_cores = 20)

# Then run for all other genomic features (on table without genes in multiple COG categories split).
combined$proMGE <- factor(combined$proMGE, levels = c('No', 'Yes'))
combined$crispr_hit <- factor(combined$crispr_hit, levels = c('No', 'Yes'))
combined$genomad <- factor(combined$genomad, levels = c('Other', 'Plasmid', 'Virus'))
combined$scaffold_contains_provirus <- factor(combined$scaffold_contains_provirus, levels = c('No', 'Yes'))

RANGER_proMGE_glmm_out <- run_glmms_per_grouping(in_tab = combined,
                                             subset_col_i = RANGER_subset_col_i,
                                             dependent_var_string = "proMGE + (1 | orderedNorm_panaroo_sp_gf_tallies) + (1 | ranger_species) + (1 | panaroo_sp_gf)",
                                             backup1_dependent_var_string = "proMGE + (1 | orderedNorm_panaroo_sp_gf_tallies) + (1 | ranger_species)",
                                             backup2_dependent_var_string = "proMGE + (1 | orderedNorm_panaroo_sp_gf_tallies)",
                                             convert_to_integer_no_yes=TRUE,
                                             num_cores = 20)

RANGER_crispr_hit_glmm_out <- run_glmms_per_grouping(in_tab = combined,
                                                subset_col_i = RANGER_subset_col_i,
                                                dependent_var_string = "crispr_hit + (1 | orderedNorm_panaroo_sp_gf_tallies) + (1 | ranger_species) + (1 | panaroo_sp_gf)",
                                                backup1_dependent_var_string = "crispr_hit + (1 | orderedNorm_panaroo_sp_gf_tallies) + (1 | ranger_species)",
                                                backup2_dependent_var_string = "crispr_hit + (1 | orderedNorm_panaroo_sp_gf_tallies)",
                                                convert_to_integer_no_yes=TRUE,
                                                num_cores = 20)

RANGER_genomad_glmm_out <- run_glmms_per_grouping(in_tab = combined,
                                                    subset_col_i = RANGER_subset_col_i,
                                                    dependent_var_string = "genomad + (1 | orderedNorm_panaroo_sp_gf_tallies) + (1 | ranger_species) + (1 | panaroo_sp_gf)",
                                                    backup1_dependent_var_string = "genomad + (1 | orderedNorm_panaroo_sp_gf_tallies) + (1 | ranger_species)",
                                                    backup2_dependent_var_string = "genomad + (1 | orderedNorm_panaroo_sp_gf_tallies)",
                                                    convert_to_integer_no_yes=TRUE,
                                                    num_cores = 20)

RANGER_scaffold_contains_provirus_glmm_out <- run_glmms_per_grouping(in_tab = combined,
                                                     subset_col_i = RANGER_subset_col_i,
                                                     dependent_var_string = "scaffold_contains_provirus + (1 | orderedNorm_panaroo_sp_gf_tallies) + (1 | ranger_species) + (1 | panaroo_sp_gf)",
                                                     backup1_dependent_var_string = "scaffold_contains_provirus + (1 | orderedNorm_panaroo_sp_gf_tallies) + (1 | ranger_species)",
                                                     backup2_dependent_var_string = "scaffold_contains_provirus + (1 | orderedNorm_panaroo_sp_gf_tallies)",
                                                     convert_to_integer_no_yes=TRUE,
                                                     num_cores = 20)

model_summaries <- list(COG_category = lapply(RANGER_COG_glmm_out, summary),
                        proMGE = lapply(RANGER_proMGE_glmm_out, summary),
                        crispr = lapply(RANGER_crispr_hit_glmm_out, summary),
                        genomad = lapply(RANGER_genomad_glmm_out, summary),
                        scaffold_w_prophage = lapply(RANGER_scaffold_contains_provirus_glmm_out, summary))

saveRDS(object = model_summaries,
        file = '/mfs/gdouglas/projects/ocean_mags/glmm_working/RANGER_glmm_model_summaries.rds')
