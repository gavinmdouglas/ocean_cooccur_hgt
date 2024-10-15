rm(list = ls(all.names = TRUE))

# Prep tables with samples split by filter cutoff.
source('~/scripts/ocean_mag_hgt/scripts/coverm/prep_coverage_tables_functions.R')

freeliv_meta <- read.table('/mfs/gdouglas/projects/ocean_mags/metadata/OceanDNA_supp_metadata/OceanDNA_metadata_water_freeliving.tsv',
                           header=TRUE, sep='\t', stringsAsFactors = FALSE)

lessfiltered_meta <- read.table('/mfs/gdouglas/projects/ocean_mags/metadata/OceanDNA_supp_metadata/OceanDNA_metadata_water_notparticledepleted.tsv',
                                header=TRUE, sep='\t', stringsAsFactors = FALSE)

coverm_set1 <- coverm_read_by_folder("/mfs/gdouglas/projects/ocean_mags/coverm/additional_output")
coverm_set2 <- coverm_read_by_folder("/mfs/gdouglas/projects/ocean_mags/coverm/additional_OceanDNA_round2")

# Confirm that all expected metagenomics samples are present.
# First, run sanity check that no samples intersect between these two sets.
if (length(intersect(unique(coverm_set1$mgs_sample), unique(coverm_set2$mgs_sample))) > 0) {
  stop('ERROR - metagenomics sample output in both sets!')
}

coverm_all <- rbind(coverm_set1, coverm_set2)
rownames(coverm_all) <- NULL

# First, confirm that all Tara ocean (DNA) samples are present.
PRJEB1787_metadata <- read.table('~/projects/ocean_mags/metadata/PRJEB1787_metadata.csv',
                                 header = TRUE, sep = ',', stringsAsFactors = FALSE)
if (length(which(! PRJEB1787_metadata$Sample %in% coverm_all$mgs_sample)) > 0) {
  stop('ERROR - some samples in PRJEB1787 missing!')
} else {
  number_PRJEB1787_present <- length(which(PRJEB1787_metadata$Sample %in% coverm_all$mgs_sample))
  message(number_PRJEB1787_present, ' of ', nrow(PRJEB1787_metadata), ' PRJEB1787 samples present in CoverM output.')
}

PRJEB9740_metadata <- read.table('~/projects/ocean_mags/metadata/PRJEB9740_metadata.csv',
                                 header = TRUE, sep = ',', stringsAsFactors = FALSE)
if (length(which(! PRJEB9740_metadata$Sample %in% coverm_all$mgs_sample)) > 0) {
  stop('ERROR - some samples in PRJEB9740 missing!')
} else {
 number_PRJEB9740_present <- length(which(PRJEB9740_metadata$Sample %in% coverm_all$mgs_sample))
 message(number_PRJEB9740_present, ' of ', nrow(PRJEB9740_metadata), ' PRJEB9740 samples present in CoverM output.')
}

# Note that the above checks are mainly sanity checks, as all but one
# Tara sample (ERS488919) is present in the OceanDNA water sample metadata anyway.
OceanDNA_water_metadata <- read.table('/mfs/gdouglas/projects/ocean_mags/metadata/OceanDNA_supp_metadata/Supp_File_S1_water_samples.tsv',
                                      header = TRUE, sep = '\t', stringsAsFactors = FALSE)
OceanDNA_water_metadata$sample_name[which(OceanDNA_water_metadata$sample_name == 'ERS492821_ERS492814')] <- 'ERS492814'
if (length(which(! OceanDNA_water_metadata$sample_name %in% coverm_all$mgs_sample)) > 0) {
  stop('ERROR - some samples in PRJEB9740 missing!')
} else {
  number_OceanDNA_present <- length(which(OceanDNA_water_metadata$sample_name %in% coverm_all$mgs_sample))
  message(number_OceanDNA_present, ' of ', nrow(OceanDNA_water_metadata), ' OceanDNA samples present in CoverM output.')
}
message('This is the number of unique CoverM samples (pre-filtering): ', length(unique(coverm_all$mgs_sample)))

# Then get table split by each subset.
coverm_freeliv <- coverm_all
coverm_freeliv <- coverm_freeliv[which(coverm_freeliv$mgs_sample %in% freeliv_meta$sample_name), ]

# Create presence tables.
coverm_freeliv_presence <- prep_presence_tab(cov_tab = coverm_freeliv,
                                             presence_cutoff = 0.30)

# Create RPKM tables, based on the same samples/genomes as in the presence tables.
coverm_freeliv_rpkm <- reshape2::dcast(data = coverm_freeliv,
                                       formula = mgs_sample ~ genome,
                                       value.var = "rpkm",
                                       fill = 0)

rownames(coverm_freeliv_rpkm) <- coverm_freeliv_rpkm$mgs_sample
coverm_freeliv_rpkm <- coverm_freeliv_rpkm[, -which(colnames(coverm_freeliv_rpkm) == 'mgs_sample')]

coverm_freeliv_rpkm <- coverm_freeliv_rpkm[rownames(coverm_freeliv_presence), colnames(coverm_freeliv_presence)]

# Also, set RPKM values for taxa with breadth below the presence cut-off to 0.
coverm_freeliv_rpkm[coverm_freeliv_presence == 0] <- 0

write_gzip_table_w_rowcol(x = coverm_freeliv_presence,
                          outfile = '~/projects/ocean_mags/networks/combined_tables/metaG_presence_freeliv.tsv.gz')

write_gzip_table_w_rowcol(x = coverm_freeliv_rpkm,
                          outfile = '~/projects/ocean_mags/networks/combined_tables/metaG_rpkm_freeliv.tsv.gz')


# And for "lessfiltered" set.
coverm_lessfiltered <- coverm_all
coverm_lessfiltered <- coverm_lessfiltered[which(coverm_lessfiltered$mgs_sample %in% lessfiltered_meta$sample_name), ]

# Create presence tables.
coverm_lessfiltered_presence <- prep_presence_tab(cov_tab = coverm_lessfiltered,
                                             presence_cutoff = 0.30)

# Create RPKM tables, based on the same samples/genomes as in the presence tables.
coverm_lessfiltered_rpkm <- reshape2::dcast(data = coverm_lessfiltered,
                                       formula = mgs_sample ~ genome,
                                       value.var = "rpkm",
                                       fill = 0)

rownames(coverm_lessfiltered_rpkm) <- coverm_lessfiltered_rpkm$mgs_sample
coverm_lessfiltered_rpkm <- coverm_lessfiltered_rpkm[, -which(colnames(coverm_lessfiltered_rpkm) == 'mgs_sample')]

coverm_lessfiltered_rpkm <- coverm_lessfiltered_rpkm[rownames(coverm_lessfiltered_presence), colnames(coverm_lessfiltered_presence)]

# Also, set RPKM values for taxa with breadth below the presence cut-off to 0.
coverm_lessfiltered_rpkm[coverm_lessfiltered_presence == 0] <- 0

write_gzip_table_w_rowcol(x = coverm_lessfiltered_presence,
                          outfile = '~/projects/ocean_mags/networks/combined_tables/metaG_presence_lessfiltered.tsv.gz')

write_gzip_table_w_rowcol(x = coverm_lessfiltered_rpkm,
                          outfile = '~/projects/ocean_mags/networks/combined_tables/metaG_rpkm_lessfiltered.tsv.gz')
