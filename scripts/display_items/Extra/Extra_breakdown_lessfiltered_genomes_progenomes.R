rm(list = ls(all.names = TRUE))

taxonomy <- read.table('/mfs/gdouglas/projects/ocean_mags/progenomes_analyses/genome_prep/progenomes_ncbi_taxonomy.tsv.gz',
                       header=TRUE, sep = '\t', stringsAsFactors = FALSE, row.names=1)

prepped <- read.table('~/projects/ocean_mags/progenomes_analyses/networks/allsamples/clusterbased/metaG_hyperg_cooccur.allsamples.combined.prepped.tsv.gz',
                      header=TRUE, sep='\t', stringsAsFactors = FALSE)

genomes_mainly_freeliving <- read.table(file = '/mfs/gdouglas/projects/ocean_hgt_zenodo/taxa_subsets/progenomes_freeliving_associated.tsv.gz',
                                        header=FALSE, stringsAsFactors = FALSE)$V1

genomes_mainly_lessfiltered <- read.table(file = '/mfs/gdouglas/projects/ocean_hgt_zenodo/taxa_subsets/progenomes_lessfiltered_associated.tsv.gz',
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

length(unique(c(prepped_lessfilt_subset$taxon_i, prepped_lessfilt_subset$taxon_j)))

length(unique(c(prepped_freeliv_subset$taxon_i, prepped_freeliv_subset$taxon_j)))

tmp <- readRDS('projects/ocean_mags/glmm_working/combined_summary.rds')
