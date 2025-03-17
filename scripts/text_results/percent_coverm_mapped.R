rm(list = ls(all.names = TRUE))

source('~/scripts/ocean_mag_hgt/scripts/coverm/prep_coverage_tables_functions.R')

coverm_set1 <- coverm_read_by_folder("/mfs/gdouglas/projects/ocean_mags/coverm/additional_output")
coverm_set2 <- coverm_read_by_folder("/mfs/gdouglas/projects/ocean_mags/coverm/additional_OceanDNA_round2")

# Confirm that all expected metagenomics samples are present.
# First, run sanity check that no samples intersect between these two sets.
if (length(intersect(unique(coverm_set1$mgs_sample), unique(coverm_set2$mgs_sample))) > 0) {
  stop('ERROR - metagenomics sample output in both sets!')
}

coverm_all <- rbind(coverm_set1, coverm_set2)
rownames(coverm_all) <- NULL

# Exclude samples that were filtered out due to having very few genomes called as present.
presence_tab <- read.table('~/projects/ocean_mags/networks/combined_tables/metaG_presence_allsamples.tsv.gz',
                           header=TRUE, sep = '\t', stringsAsFactors = FALSE, row.names=1)

coverm_all <- coverm_all[which(coverm_all$mgs_sample %in% rownames(presence_tab)), ]

length(unique(coverm_all$mgs_sample))
# 1862

coverm_mapped <- coverm_all[-which(is.na(coverm_all$readcount)), ]

# Mean and SD number of reads mapped per sample.
sum_mapped_counts <- aggregate(readcount ~ mgs_sample, FUN = sum, data = coverm_mapped)
mean(sum_mapped_counts$readcount)
# 26592809
sd(sum_mapped_counts$readcount)
# 36857948

unmapped <- coverm_all[which(coverm_all$genome == 'unmapped'), ]

mapped_percent <- 100.0 - unmapped$relabun

mean(mapped_percent)
# 23.88197

sd(mapped_percent)
# 10.1516

