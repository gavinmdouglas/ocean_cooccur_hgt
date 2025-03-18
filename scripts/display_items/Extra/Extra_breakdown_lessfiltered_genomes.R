rm(list = ls(all.names = TRUE))

taxonomy <- read.table('/mfs/gdouglas/projects/ocean_hgt_zenodo/mapfiles/MAG_taxa_breakdown.tsv.gz',
                       header=TRUE, sep = '\t', stringsAsFactors = FALSE, row.names=2)

prepped <- read.table('~/projects/ocean_mags/networks/allsamples/clusterbased/metaG_hyperg_cooccur.allsamples.combined.prepped.tsv.gz',
                      header=TRUE, sep='\t', stringsAsFactors = FALSE)

genomes_mainly_freeliving <- read.table(file = '/mfs/gdouglas/projects/ocean_hgt_zenodo/taxa_subsets/Tara_genomes_freeliving_associated.tsv.gz',
                                        header=FALSE, stringsAsFactors = FALSE)$V1

genomes_mainly_lessfiltered <- read.table(file = '/mfs/gdouglas/projects/ocean_hgt_zenodo/taxa_subsets/Tara_genomes_lessfiltered_associated.tsv.gz',
                                          header=FALSE, stringsAsFactors = FALSE)$V1

taxon_i_grouping <- rep('Other', nrow(prepped))
taxon_j_grouping <- rep('Other', nrow(prepped))
taxon_i_grouping[which(prepped$taxon_i %in% genomes_mainly_freeliving)] <- 'freeliving'
taxon_i_grouping[which(prepped$taxon_i %in% genomes_mainly_lessfiltered)] <- 'lessfiltered'
taxon_j_grouping[which(prepped$taxon_j %in% genomes_mainly_freeliving)] <- 'freeliving'
taxon_j_grouping[which(prepped$taxon_j %in% genomes_mainly_lessfiltered)] <- 'lessfiltered'

prepped$filter_group_match <- 'Mixed'
prepped[which(taxon_i_grouping == 'freeliving' & taxon_j_grouping == 'freeliving'), 'filter_group_match'] <- 'Free-living'
prepped[which(taxon_i_grouping == 'lessfiltered' & taxon_j_grouping == 'lessfiltered'), 'filter_group_match'] <- 'Less-filtered'


prepped_freeliv_subset <- prepped[which(prepped$filter_group_match == 'Free-living'), ]
prepped_lessfilt_subset <- prepped[which(prepped$filter_group_match == 'Less-filtered'), ]

prepped_freeliv_hgt_subset <- prepped_freeliv_subset[which(prepped_freeliv_subset$hgt == 1), ]
prepped_lessfilt_hgt_subset <- prepped_lessfilt_subset[which(prepped_lessfilt_subset$hgt == 1), ]

length(unique(c(prepped_freeliv_hgt_subset$taxon_i, prepped_freeliv_hgt_subset$taxon_j)))
length(unique(c(prepped_lessfilt_hgt_subset$taxon_i, prepped_lessfilt_hgt_subset$taxon_j)))

lessfilt_hgt_taxa <- unique(c(prepped_lessfilt_hgt_subset$taxon_i, prepped_lessfilt_hgt_subset$taxon_j))
lessfilt_hgt_taxonomy <- taxonomy[lessfilt_hgt_taxa, ]

write.table(x = lessfilt_hgt_taxonomy, file = '~/tmp/lessfilt_hgt_taxonomy.tsv', col.names = NA, row.names = TRUE, sep = '\t', quote = FALSE)

write.table(x = prepped_lessfilt_hgt_subset, file = '~/tmp/lessfilt_hgt_prepped_subset.tsv', col.names = TRUE, row.names = FALSE, sep = '\t', quote = FALSE)
