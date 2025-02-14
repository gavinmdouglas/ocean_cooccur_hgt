rm(list = ls(all.names = TRUE))

# Prep separate tables for samples in Tara and GEOTRACES datasets specifically.
source('~/scripts/ocean_mag_hgt/scripts/coverm/prep_coverage_tables_functions.R')

OceanDNA_water_metadata <- read.table('/mfs/gdouglas/projects/ocean_mags/metadata/OceanDNA_supp_metadata/Supp_File_S1_water_samples.tsv',
                                      header = TRUE, sep = '\t', stringsAsFactors = FALSE)

geotraces_meta <- OceanDNA_water_metadata[which(OceanDNA_water_metadata$division == 'GEOTRACES'), ]
tara_meta <- OceanDNA_water_metadata[grep('^Tara', OceanDNA_water_metadata$division), ]

tara_meta[which(tara_meta$sample_name == 'ERS492821_ERS492814'), 'sample_name'] <- 'ERS492814'

coverm_set1 <- coverm_read_by_folder("/mfs/gdouglas/projects/ocean_mags/coverm/additional_output")
coverm_set2 <- coverm_read_by_folder("/mfs/gdouglas/projects/ocean_mags/coverm/additional_OceanDNA_round2")

# Confirm that all expected metagenomics samples are present.
# First, run sanity check that no samples intersect between these two sets.
if (length(intersect(unique(coverm_set1$mgs_sample), unique(coverm_set2$mgs_sample))) > 0) {
  stop('ERROR - metagenomics sample output in both sets!')
}

coverm_all <- rbind(coverm_set1, coverm_set2)
rownames(coverm_all) <- NULL

coverm_geotraces <- coverm_all[which(coverm_all$mgs_sample %in% geotraces_meta$sample_name), ]
coverm_tara <- coverm_all[which(coverm_all$mgs_sample %in% tara_meta$sample_name), ]

# Create presence tables.
coverm_geotraces_presence <- prep_presence_tab(cov_tab = coverm_geotraces,
                                               presence_cutoff = 0.30)

# Create RPKM tables, based on the same samples/genomes as in the presence tables.
coverm_geotraces_rpkm <- reshape2::dcast(data = coverm_geotraces,
                                       formula = mgs_sample ~ genome,
                                       value.var = "rpkm",
                                       fill = 0)

rownames(coverm_geotraces_rpkm) <- coverm_geotraces_rpkm$mgs_sample
coverm_geotraces_rpkm <- coverm_geotraces_rpkm[, -which(colnames(coverm_geotraces_rpkm) == 'mgs_sample')]

coverm_geotraces_rpkm <- coverm_geotraces_rpkm[rownames(coverm_geotraces_presence), colnames(coverm_geotraces_presence)]

# Also, set RPKM values for taxa with breadth below the presence cut-off to 0.
coverm_geotraces_rpkm[coverm_geotraces_presence == 0] <- 0

write_gzip_table_w_rowcol(x = coverm_geotraces_presence,
                          outfile = '~/projects/ocean_mags/networks/combined_tables/metaG_presence_geotraces.tsv.gz')

write_gzip_table_w_rowcol(x = coverm_geotraces_rpkm,
                          outfile = '~/projects/ocean_mags/networks/combined_tables/metaG_rpkm_geotraces.tsv.gz')


# And for "tara" set.
coverm_tara <- coverm_all
coverm_tara <- coverm_tara[which(coverm_tara$mgs_sample %in% tara_meta$sample_name), ]

# Create presence tables.
coverm_tara_presence <- prep_presence_tab(cov_tab = coverm_tara,
                                             presence_cutoff = 0.30)

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
                          outfile = '~/projects/ocean_mags/networks/combined_tables/metaG_presence_tara.tsv.gz')

write_gzip_table_w_rowcol(x = coverm_tara_rpkm,
                          outfile = '~/projects/ocean_mags/networks/combined_tables/metaG_rpkm_tara.tsv.gz')
