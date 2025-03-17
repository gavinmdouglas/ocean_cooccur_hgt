# Run COG category enrichments for 95-<99% and >=99% putatively transferred HGT hits at each taxonomic level.
# But this time, only consider COG-annotated clusters.

rm(list = ls(all.names = TRUE))

source('~/scripts/ocean_mag_hgt/scripts/func_enrich/enrich_function.R')

COG_category_descrip <- read.table("/mfs/gdouglas/db/COG_definitions/COG_category_descrip.tsv",
                                   header = FALSE, sep = "\t", row.names = 1)

cluster_annot <- read.table(file = "/mfs/gdouglas/projects/ocean_mags/Sunagawa_dataset/genomes-cdhit-cluster-cog-info.tsv.gz",
                            sep = "\t", row.names = 1, stringsAsFactors = FALSE, header = TRUE)

gene_to_cluster <- read.table("/mfs/gdouglas/projects/ocean_mags/clusters/gene_to_cluster.tsv.gz",
                              header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = 1)


best_hits <- read.table(file = '/mfs/gdouglas/projects/ocean_mags/clusters/all_best_hits.tsv.gz', header = TRUE,
                        sep = '\t', stringsAsFactors = FALSE)

best_hits$gene1_cluster <- gene_to_cluster[best_hits$gene1, 'cluster_id']
best_hits$gene2_cluster <- gene_to_cluster[best_hits$gene2, 'cluster_id']

# Remove unannot clusters.
cluster_annot <- cluster_annot[which(cluster_annot$majority_rule_COG_category != "-"), ]
cluster_annot_mapping <- list()
background_all_genomes <- as.character()
for (COG_category in rownames(COG_category_descrip)) {
  cluster_annot_mapping[[COG_category]] <- rownames(cluster_annot)[grep(COG_category, cluster_annot$majority_rule_COG_category)]
  background_all_genomes <- unique(c(background_all_genomes, rownames(cluster_annot)[grep(COG_category, cluster_annot$majority_rule_COG_category)]))
}

best_hits <- best_hits[which(best_hits$gene1_cluster %in% rownames(cluster_annot)), ]
best_hits <- best_hits[which(best_hits$gene2_cluster %in% rownames(cluster_annot)), ]

to_ignore <- c("A", "B", "Y", "Z")

all_output_raw <- list()

identities <- c('95', '99', 'both')
taxa_levels <- c('Genus', 'Family', 'Order', 'Class', 'Phylum', 'Domain')

best_hits_by_tax <- list()
best_hits_by_tax[['95']] <- list()
best_hits_by_tax[['99']] <- list()
best_hits_by_tax[['both']] <- list()

for (tax_level in taxa_levels) {
  best_hit_subset <- best_hits[which(best_hits$highest_tax_diff == tax_level), ]
  best_hits_by_tax[['both']][[tax_level]] <- unique(c(best_hit_subset$gene1_cluster, best_hit_subset$gene2_cluster))

  best_hit_subset_95 <- best_hit_subset[which(best_hit_subset$identity < 99 & best_hit_subset$identity >= 95), ]
  if (nrow(best_hit_subset_95) > 0) {
    best_hits_by_tax[['95']][[tax_level]] <- unique(c(best_hit_subset_95$gene1_cluster, best_hit_subset_95$gene2_cluster))
  } else {
    best_hits_by_tax[['95']][[tax_level]] <- character()
  }

  best_hit_subset_99 <- best_hit_subset[which(best_hit_subset$identity >= 99), ]
  if (nrow(best_hit_subset_99) > 0) {
    best_hits_by_tax[['99']][[tax_level]] <- unique(c(best_hit_subset_99$gene1_cluster, best_hit_subset_99$gene2_cluster))
  } else {
    best_hits_by_tax[['99']][[tax_level]] <- character()
  }
}

combo_num <- 1
for (identity in identities) {
  for (taxa_level in taxa_levels) {
    raw_level_out <- identify_enriched_categories(genes = best_hits_by_tax[[identity]][[taxa_level]],
                                                  background = background_all_genomes[which(! background_all_genomes %in% best_hits_by_tax[[identity]][[taxa_level]])],
                                                  gene_to_category_map = cluster_annot_mapping,
                                                  to_ignore = to_ignore)
    raw_level_out$identity <- identity
    raw_level_out$taxon <- taxa_level
    all_output_raw[[combo_num]] <- raw_level_out
    combo_num <- combo_num + 1
  }
}

all_output <- do.call(rbind, all_output_raw)
all_output$descrip <- COG_category_descrip[all_output$category, 1]

all_output$fdr <- p.adjust(all_output$p, 'BH')

write.table(x = all_output, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE,
            file = "/mfs/gdouglas/projects/ocean_hgt_zenodo/putative_hgt/cluster/clusterbased_COG_category_enrichment_no.unannot.tsv")
