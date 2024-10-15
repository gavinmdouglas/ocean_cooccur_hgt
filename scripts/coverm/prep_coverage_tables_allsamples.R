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

# Create presence tables.
dna_presence <- prep_presence_tab(cov_tab = coverm_all,
                                  presence_cutoff = 0.30)

# Create RPKM tables, based on the same samples/genomes as in the presence tables.
dna_rpkm <- reshape2::dcast(data = coverm_all,
                            formula = mgs_sample ~ genome,
                            value.var = "rpkm",
                            fill = 0)

rownames(dna_rpkm) <- dna_rpkm$mgs_sample
dna_rpkm <- dna_rpkm[, -which(colnames(dna_rpkm) == 'mgs_sample')]

dna_rpkm <- dna_rpkm[rownames(dna_presence), colnames(dna_presence)]

# Also, set RPKM values for taxa with breadth below the presence cut-off to 0.
dna_rpkm[dna_presence == 0] <- 0

# Check whether ERS488919 (non-OceanDNA sample) filtered out due to prevalence cut-offs anyway.
if ("ERS488919" %in% rownames(dna_presence)) {
  message('Extra Tara ocean sample made it through filters...')
}

write_gzip_table_w_rowcol(x = dna_presence,
                          outfile = '~/projects/ocean_mags/networks/combined_tables/metaG_presence_allsamples.tsv.gz')

write_gzip_table_w_rowcol(x = dna_rpkm,
                          outfile = '~/projects/ocean_mags/networks/combined_tables/metaG_rpkm_allsamples.tsv.gz')
