rm(list = ls(all.names = TRUE))

source('~/scripts/ocean_cooccur_hgt/scripts/coverm/prep_coverage_tables_functions.R')

coverm_all <- coverm_read_by_folder("/mfs/gdouglas/projects/ocean_mags/additional_protist_fraction/coverm/output")

# First, confirm that all Tara ocean (DNA) samples are present.
PRJEB4352_metadata <- read.table('~/projects/ocean_mags/additional_protist_fraction/metadata/PRJEB4352_ENA_metadata_filt_focal.txt',
                                 header = TRUE, sep = '\t', stringsAsFactors = FALSE)
if (length(which(! PRJEB4352_metadata$sample_accession %in% coverm_all$mgs_sample)) > 0) {
  stop('ERROR - some samples in PRJEB4352 missing!')
}

PRJEB9691_metadata <- read.table('~/projects/ocean_mags/additional_protist_fraction/metadata/PRJEB9691_ENA_metadata_filt_focal.txt',
                                 header = TRUE, sep = '\t', stringsAsFactors = FALSE)
if (length(which(! PRJEB9691_metadata$sample_accession %in% coverm_all$mgs_sample)) > 0) {
  stop('ERROR - some samples in PRJEB9691 missing!')
}

dna_presence <- prep_presence_tab(cov_tab = coverm_all,
                                  presence_cutoff = 0.30,
                                  min_samples_w_genome = 2,
                                  min_genomes_per_sample = 1)

protist_assoc_taxa <- colnames(dna_presence)

write.table(x = protist_assoc_taxa,
            file = '~/tmp/protist_assoc_taxa.txt',
            sep = '\t', col.names = FALSE, row.names = FALSE, quote = FALSE)
