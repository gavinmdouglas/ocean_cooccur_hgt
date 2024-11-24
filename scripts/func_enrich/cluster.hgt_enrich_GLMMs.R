rm(list = ls(all.names = TRUE))

source('/mfs/gdouglas/scripts/ocean_mag_hgt/scripts/functions.R')

library(glmmTMB)

combined <- read.table('/mfs/gdouglas/projects/ocean_mags/summary_files/gene_info_and_hgt.tsv.gz',
                       header = TRUE, sep = '\t', stringsAsFactors = FALSE)

# First run test on COG categories.
combined_COG <- combined
combined_COG <- combined_COG[which(combined_COG$COG_category != '-'), ]
combined_COG <- split_multi_category_rows(in_df = combined_COG, category_col = 'COG_category', num_cores = 100)
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

Clusterbased_subset_col_i <- grep("Clusterbased", colnames(combined_COG))

Clusterbased_COG_glmm_out <- run_glmms_per_grouping(in_tab = combined_COG,
                               subset_col_i = Clusterbased_subset_col_i,
                               dependent_var_string = "COG_category + (1 | gene_rep) + (1 | genome)",
                               backup1_dependent_var_string = "COG_category + (1 | genome)",
                               backup2_dependent_var_string = "COG_category",
                               convert_to_integer_no_yes=TRUE,
                               num_cores = 100)


# Then run for all other genomic features (on table without genes in multiple COG categories split).
combined$proMGE <- factor(combined$proMGE, levels = c('No', 'Yes'))
combined$crispr_hit <- factor(combined$crispr_hit, levels = c('No', 'Yes'))
combined$genomad <- factor(combined$genomad, levels = c('Other', 'Plasmid', 'Virus'))
combined$scaffold_contains_provirus <- factor(combined$scaffold_contains_provirus, levels = c('No', 'Yes'))

Clusterbased_proMGE_glmm_out <- run_glmms_per_grouping(in_tab = combined,
                                             subset_col_i = Clusterbased_subset_col_i,
                                             dependent_var_string = "proMGE + (1 | gene_rep) + (1 | genome)",
                                             backup1_dependent_var_string = "proMGE + (1 | genome)",
                                             backup2_dependent_var_string = "proMGE",
                                             convert_to_integer_no_yes=TRUE,
                                             num_cores = 100)

Clusterbased_crispr_hit_glmm_out <- run_glmms_per_grouping(in_tab = combined,
                                                subset_col_i = Clusterbased_subset_col_i,
                                                dependent_var_string = "crispr_hit + (1 | gene_rep) + (1 | genome)",
                                                backup1_dependent_var_string = "crispr_hit + (1 | genome)",
                                                backup2_dependent_var_string = "crispr_hit",
                                                convert_to_integer_no_yes=TRUE,
                                                num_cores = 100)

Clusterbased_genomad_glmm_out <- run_glmms_per_grouping(in_tab = combined,
                                                    subset_col_i = Clusterbased_subset_col_i,
                                                    dependent_var_string = "genomad + (1 | gene_rep) + (1 | genome)",
                                                    backup1_dependent_var_string = "genomad + (1 | genome)",
                                                    backup2_dependent_var_string = "genomad",
                                                    convert_to_integer_no_yes=TRUE,
                                                    num_cores = 100)

Clusterbased_scaffold_contains_provirus_glmm_out <- run_glmms_per_grouping(in_tab = combined,
                                                     subset_col_i = Clusterbased_subset_col_i,
                                                     dependent_var_string = "scaffold_contains_provirus + (1 | gene_rep) + (1 | genome)",
                                                     backup1_dependent_var_string = "scaffold_contains_provirus + (1 | genome)",
                                                     backup2_dependent_var_string = "scaffold_contains_provirus",
                                                     convert_to_integer_no_yes=TRUE,
                                                     num_cores = 100)

model_summaries <- list(COG_category = lapply(Clusterbased_COG_glmm_out, summary),
                        proMGE = lapply(Clusterbased_proMGE_glmm_out, summary),
                        crispr = lapply(Clusterbased_crispr_hit_glmm_out, summary),
                        genomad = lapply(Clusterbased_genomad_glmm_out, summary),
                        scaffold_w_prophage = lapply(Clusterbased_scaffold_contains_provirus_glmm_out, summary))

saveRDS(object = model_summaries,
        file = '/mfs/gdouglas/projects/ocean_mags/glmm_working/Clusterbased_glmm_model_summaries.rds')
