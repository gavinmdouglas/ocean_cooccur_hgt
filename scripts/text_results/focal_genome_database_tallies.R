rm(list = ls(all.names = TRUE))

taxonomy <- read.table('/mfs/gdouglas/projects/ocean_hgt_zenodo/mapfiles/MAG_taxa_breakdown.tsv.gz',
                       header=TRUE, sep = '\t', stringsAsFactors = FALSE, row.names=2)

# Total
nrow(taxonomy)

# Number of MAGs
mag_ids <- grep("_METAG_", taxonomy$MAG, value = TRUE)
length(mag_ids)
if ((length(grep("_REFG_", mag_ids)) > 0) | (length(grep("_SAGS_", mag_ids)) > 0)) {
  print('Other hits too!')
}

# Number of reference genomes
ref_genome_ids <- grep("_REFG_", taxonomy$MAG, value = TRUE)
length(ref_genome_ids)
if ((length(grep("_METAG_", ref_genome_ids)) > 0) | (length(grep("_SAGS_", ref_genome_ids)) > 0)) {
  print('Other hits too!')
}

# Number of SAGS
sag_ids <- grep("_SAGS_", taxonomy$MAG, value = TRUE)
length(sag_ids)
if ((length(grep("_METAG_", sag_ids)) > 0) | (length(grep("_METAG_", sag_ids)) > 0)) {
  print('Other hits too!')
}