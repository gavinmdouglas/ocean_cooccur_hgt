rm(list = ls(all.names = TRUE))

combined_summary <- read.table('/mfs/gdouglas/projects/ocean_hgt_zenodo/hgt_prev_analyses/Tara_fraction_sample_matched_HGT_prev.tsv.gz',
                               header=TRUE, sep = '\t', stringsAsFactors = FALSE)

combined_summary <- combined_summary[which(combined_summary$lower_total_genomes >= 10 & combined_summary$upper_total_genomes >= 10), ]

# Also decided to remove all "prop_mult" panels, it we decided this was too much information,
# and unclear how to interpret (at the moment at least!)

col_to_rm <- grep("filt_group|_sample|total_genomes|mean_genes|_prop_mult_", colnames(combined_summary), value = TRUE)
combined_summary <- combined_summary[, which(! colnames(combined_summary) %in% col_to_rm)]

mean(combined_summary$upper_prop_hgt / combined_summary$lower_prop_hgt)
sd(combined_summary$upper_prop_hgt / combined_summary$lower_prop_hgt)
