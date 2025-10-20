rm(list = ls(all.names = TRUE))

# Prep coverage tables as before, but include the protist-fraction samples as well.
# Create the full table for all samples (with protist-fraction), as well as Tara ocean samples only.

source('~/scripts/ocean_cooccur_hgt/scripts/coverm/prep_coverage_tables_functions.R')

coverm_protist <- coverm_read_by_folder("/mfs/gdouglas/projects/ocean_mags/additional_protist_fraction/coverm/output")

coverm_protist$subfolder <- "protist_fraction"

# First, confirm that all Tara ocean (DNA) samples are present.
PRJEB4352_metadata <- read.table('~/projects/ocean_mags/additional_protist_fraction/metadata/PRJEB4352_ENA_metadata_filt_focal.txt',
                                 header = TRUE, sep = '\t', stringsAsFactors = FALSE)
if (length(which(! PRJEB4352_metadata$sample_accession %in% coverm_protist$mgs_sample)) > 0) {
  stop('ERROR - some samples in PRJEB4352 missing!')
}

PRJEB9691_metadata <- read.table('~/projects/ocean_mags/additional_protist_fraction/metadata/PRJEB9691_ENA_metadata_filt_focal.txt',
                                 header = TRUE, sep = '\t', stringsAsFactors = FALSE)
if (length(which(! PRJEB9691_metadata$sample_accession %in% coverm_protist$mgs_sample)) > 0) {
  stop('ERROR - some samples in PRJEB9691 missing!')
}

coverm_set1 <- coverm_read_by_folder("/mfs/gdouglas/projects/ocean_mags/coverm/additional_output")
coverm_set2 <- coverm_read_by_folder("/mfs/gdouglas/projects/ocean_mags/coverm/additional_OceanDNA_round2")

# Confirm that all expected metagenomics samples are present.
# First, run sanity check that no samples intersect between these two sets.
if (length(intersect(unique(coverm_set1$mgs_sample), unique(coverm_set2$mgs_sample))) > 0) {
  stop('ERROR - metagenomics sample output in both sets!')
}

coverm_orig <- rbind(coverm_set1, coverm_set2)
rownames(coverm_orig) <- NULL

# Get subset of Tara only.
OceanDNA_water_metadata <- read.table('/mfs/gdouglas/projects/ocean_mags/metadata/OceanDNA_supp_metadata/Supp_File_S1_water_samples.tsv',
                                      header = TRUE, sep = '\t', stringsAsFactors = FALSE)
tara_meta <- OceanDNA_water_metadata[grep('^Tara', OceanDNA_water_metadata$division), ]
tara_meta[which(tara_meta$sample_name == 'ERS492821_ERS492814'), 'sample_name'] <- 'ERS492814'

coverm_tara <- coverm_orig[which(coverm_orig$mgs_sample %in% tara_meta$sample_name), ]
coverm_tara <- rbind(coverm_tara, coverm_protist)

coverm_tara_presence <- prep_presence_tab(cov_tab = coverm_tara,
                                          presence_cutoff = 0.30)

coverm_tara_presence <- coverm_tara_presence[which(rowSums(coverm_tara_presence) > 0), which(colSums(coverm_tara_presence) > 0)]

# Create RPKM tables, based on the same samples/genomes as in the presence tables.
coverm_tara_rpkm <- reshape2::dcast(data = coverm_tara,
                                         formula = mgs_sample ~ genome,
                                         value.var = "rpkm",
                                         fill = 0)

rownames(coverm_tara_rpkm) <- coverm_tara_rpkm$mgs_sample
coverm_tara_rpkm <- coverm_tara_rpkm[, -which(colnames(coverm_tara_rpkm) == 'mgs_sample')]

coverm_tara_rpkm <- coverm_tara_rpkm[rownames(coverm_tara_presence), colnames(coverm_tara_presence)]

# Also, set RPKM values for taxa with breadth below the presence cut-off to 0.
coverm_tara_rpkm[coverm_tara_presence == 0] <- 0

write_gzip_table_w_rowcol(x = coverm_tara_presence,
                          outfile = '~/projects/ocean_mags/networks/combined_tables/metaG_presence_tara_w_protist_fraction.tsv.gz')

write_gzip_table_w_rowcol(x = coverm_tara_rpkm,
                          outfile = '~/projects/ocean_mags/networks/combined_tables/metaG_rpkm_tara_w_protist_fraction.tsv.gz')






coverm_all <- rbind(coverm_orig, coverm_protist)

coverm_all_presence <- prep_presence_tab(cov_tab = coverm_all,
                                         presence_cutoff = 0.30)

coverm_all_presence <- coverm_all_presence[which(rowSums(coverm_all_presence) > 0), which(colSums(coverm_all_presence) > 0)]

# Create RPKM tables, based on the same samples/genomes as in the presence tables.
coverm_all_rpkm <- reshape2::dcast(data = coverm_all,
                                   formula = mgs_sample ~ genome,
                                   value.var = "rpkm",
                                   fill = 0)

rownames(coverm_all_rpkm) <- coverm_all_rpkm$mgs_sample
coverm_all_rpkm <- coverm_all_rpkm[, -which(colnames(coverm_all_rpkm) == 'mgs_sample')]

coverm_all_rpkm <- coverm_all_rpkm[rownames(coverm_all_presence), colnames(coverm_all_presence)]

# Also, set RPKM values for taxa with breadth below the presence cut-off to 0.
coverm_all_rpkm[coverm_all_presence == 0] <- 0

write_gzip_table_w_rowcol(x = coverm_all_presence,
                          outfile = '~/projects/ocean_mags/networks/combined_tables/metaG_presence_all_w_protist_fraction.tsv.gz')

write_gzip_table_w_rowcol(x = coverm_all_rpkm,
                          outfile = '~/projects/ocean_mags/networks/combined_tables/metaG_rpkm_all_w_protist_fraction.tsv.gz')
