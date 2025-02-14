rm(list = ls(all.names = TRUE))

# Prep separate tables for samples in Tara and GEOTRACES datasets specifically.
source('~/scripts/ocean_mag_hgt/scripts/coverm/prep_coverage_tables_functions.R')

OceanDNA_water_metadata <- read.table('/mfs/gdouglas/projects/ocean_mags/metadata/OceanDNA_supp_metadata/Supp_File_S1_water_samples.tsv',
                                      header = TRUE, sep = '\t', stringsAsFactors = FALSE)

geotraces_meta <- OceanDNA_water_metadata[which(OceanDNA_water_metadata$division == 'GEOTRACES'), ]
tara_meta <- OceanDNA_water_metadata[grep('^Tara', OceanDNA_water_metadata$division), ]

tara_meta[which(tara_meta$sample_name == 'ERS492821_ERS492814'), 'sample_name'] <- 'ERS492814'

coverm_all <- coverm_read_by_folder("/mfs/gdouglas/projects/ocean_mags/progenomes_analyses/coverm/output/", exp_rows = 11061)
rownames(coverm_all) <- NULL

coverm_geotraces <- coverm_all[which(coverm_all$mgs_sample %in% geotraces_meta$sample_name), ]
coverm_tara <- coverm_all[which(coverm_all$mgs_sample %in% tara_meta$sample_name), ]

genomes_to_keep <- read.table('/mfs/gdouglas/projects/ocean_mags/progenomes_analyses/genome_prep/clearly_non_mags_genome_ids.txt',
                              header = FALSE, stringsAsFactors = FALSE)$V1

coverm_geotraces <- coverm_geotraces[which(coverm_geotraces$genome %in% genomes_to_keep), ]
coverm_tara <- coverm_tara[which(coverm_tara$genome %in% genomes_to_keep), ]

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
                          outfile = '~/projects/ocean_mags/progenomes_analyses/networks/combined_tables/metaG_presence_geotraces.tsv.gz')

write_gzip_table_w_rowcol(x = coverm_geotraces_rpkm,
                          outfile = '~/projects/ocean_mags/progenomes_analyses/networks/combined_tables/metaG_rpkm_geotraces.tsv.gz')

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
                          outfile = '~/projects/ocean_mags/progenomes_analyses/networks/combined_tables/metaG_presence_tara.tsv.gz')

write_gzip_table_w_rowcol(x = coverm_tara_rpkm,
                          outfile = '~/projects/ocean_mags/progenomes_analyses/networks/combined_tables/metaG_rpkm_tara.tsv.gz')
