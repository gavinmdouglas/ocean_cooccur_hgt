rm(list = ls(all.names = TRUE))

# Sanity check that HGT genes in genomes of larger filter-cutoff grouping are not obvious full of eukaryotic genes.
# One way to get at this is to see if COG categories associated with eukaryotes are more common in them.

best_hits <- read.table(file = '/mfs/gdouglas/projects/ocean_mags/clusters/all_best_hits_w_cluster_COG.tsv', header = TRUE,
                        sep = '\t', stringsAsFactors = FALSE)

lessfilt <- read.table('~/projects/ocean_hgt_zenodo/taxa_subsets/Tara_genomes_lessfiltered_associated.tsv.gz',
                       stringsAsFactors = FALSE, header=FALSE)$V1

freeliv <- read.table('~/projects/ocean_hgt_zenodo/taxa_subsets/Tara_genomes_freeliving_associated.tsv.gz',
                       stringsAsFactors = FALSE, header=FALSE)$V1

taxonomy <- read.table('/mfs/gdouglas/projects/ocean_hgt_zenodo/mapfiles/MAG_taxa_breakdown.tsv.gz',
                        header=TRUE, sep = '\t', stringsAsFactors = FALSE, row.names=1)

best_hits$taxon1 <- taxonomy[best_hits$gene1_genome, 'Taxon_ID']
best_hits$taxon2 <- taxonomy[best_hits$gene2_genome, 'Taxon_ID']

best_hits$grouping <- 'Mixed'
best_hits[which(best_hits$taxon1 %in% lessfilt & best_hits$taxon2 %in% lessfilt), 'grouping'] <- 'lessfilt'
best_hits[which(best_hits$taxon1 %in% freeliv & best_hits$taxon2 %in% freeliv), 'grouping'] <- 'freeliv'

best_hits_lessfilt <- best_hits[which(best_hits$grouping == 'lessfilt'), ]
best_hits_lessfilt <- best_hits[which(best_hits$grouping == 'freeliv'), ]
best_hits_lessfilt <- best_hits[which(best_hits$grouping == 'Mixed'), ]

#taxonomy[unique(c(best_hits_lessfilt$gene1_genome, best_hits_lessfilt$gene2_genome)), ]
