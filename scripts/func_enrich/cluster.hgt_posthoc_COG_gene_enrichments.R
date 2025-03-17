rm(list = ls(all.names = TRUE))

source('~/scripts/ocean_mag_hgt/scripts/func_enrich/enrich_function.R')

# Run enrichments for specific COG gene families in categories of interest (for specific identity / taxonomic levels).
# Rather than testing all possibilities, idea here is to dig into significant categories that showed inconsistent pattern
# across different levels, to see if we can make sense of what could be underlying these differences.

cog_descrip <- read.table('~/db/COG_definitions/cog-20.def.tab',
                          header=FALSE, sep = '\t', stringsAsFactors = FALSE, row.names=1,
                          quote='', comment.char = '')

cog_gene_info <- readRDS('~/db/COG_definitions/cog-20.to_category_collapse.rds')

categories_of_interest <- c('J', 'E', 'F', 'I', 'L', 'D', 'C')
category_results <- read.table('/mfs/gdouglas/projects/ocean_hgt_zenodo/putative_hgt/cluster/clusterbased_COG_category_enrichment_no.unannot.tsv.gz',
                               sep = '\t', stringsAsFactors = FALSE, header=TRUE)
subsets_of_interest <- category_results[which(category_results$identity == '99'), ]
subsets_of_interest <- subsets_of_interest[which(subsets_of_interest$fdr < 0.05), ]
subsets_of_interest <- subsets_of_interest[which(subsets_of_interest$category %in% categories_of_interest), ]

cluster_annot <- read.table(file = "/mfs/gdouglas/projects/ocean_mags/Sunagawa_dataset/genomes-cdhit-cluster-cog-info.tsv.gz",
                            sep = "\t", row.names = 1, stringsAsFactors = FALSE, header = TRUE)

gene_to_cluster <- read.table("/mfs/gdouglas/projects/ocean_mags/clusters/gene_to_cluster.tsv.gz",
                              header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = 1)

best_hits <- read.table(file = '/mfs/gdouglas/projects/ocean_mags/clusters/all_best_hits.tsv.gz', header = TRUE,
                        sep = '\t', stringsAsFactors = FALSE)

best_hits$gene1_cluster <- gene_to_cluster[best_hits$gene1, 'cluster_id']
best_hits$gene2_cluster <- gene_to_cluster[best_hits$gene2, 'cluster_id']

# Ignore clusters not in categories of interest.
cluster_annot <- cluster_annot[which(! cluster_annot$majority_rule_COG_category %in% categories_of_interest), ]

cluster_annot_mapping <- list()
background_all_genomes <- list()

for (COG_category in categories_of_interest) {
  cluster_annot_mapping[[COG_category]] <- list()
  background_all_genomes[[COG_category]] <- character()
  cog_genes <- cog_gene_info[grep(COG_category, cog_gene_info$category), 'COG']

  for (cog_gene in cog_genes) {
    cog_gene_pattern <- paste0(cog_gene, "(,|$)")
    cluster_cog_gene_matches <- rownames(cluster_annot)[grep(cog_gene_pattern, cluster_annot$majority_rule_COG)]
    background_all_genomes[[COG_category]] <- unique(c(background_all_genomes[[COG_category]], cluster_cog_gene_matches))
    if (length(cluster_cog_gene_matches) >= 10) {
      cluster_annot_mapping[[COG_category]][[cog_gene]] <- cluster_cog_gene_matches
    }
  }
}

best_hits <- best_hits[which(best_hits$gene1_cluster %in% rownames(cluster_annot)), ]
best_hits <- best_hits[which(best_hits$gene2_cluster %in% rownames(cluster_annot)), ]

best_hits <- best_hits[which(best_hits$identity >= 99), ]

taxa_levels <- c('Genus', 'Family', 'Order', 'Class', 'Phylum', 'Domain')
best_hits_by_tax_and_category <- list()

for (category in categories_of_interest) {
  best_hits_by_tax_and_category[[category]] <- list()
  for (tax_level in taxa_levels) {
    best_hit_subset <- best_hits[which(best_hits$highest_tax_diff == tax_level), ]
    best_hit_subset <- best_hit_subset[which(best_hit_subset$gene1_cluster %in% background_all_genomes[[category]]), ]
    best_hit_subset <- best_hit_subset[which(best_hit_subset$gene2_cluster %in% background_all_genomes[[category]]), ]
    best_hits_by_tax_and_category[[category]][[tax_level]] <- unique(c(best_hit_subset$gene1_cluster, best_hit_subset$gene2_cluster))
  }
}

all_output_raw <- list()
combo_num <- 1
for (category in categories_of_interest) {

  for (taxa_level in taxa_levels) {

    focal_set <- best_hits_by_tax_and_category[[category]][[taxa_level]]

    if (length(focal_set) < 10) {
      next
    }

    background_set <- background_all_genomes[[category]][which(! background_all_genomes[[category]] %in% focal_set)]

    raw_level_out <- identify_enriched_categories(genes = best_hits_by_tax_and_category[[category]][[taxa_level]],
                                                  background = background_set,
                                                  gene_to_category_map = cluster_annot_mapping[[category]])
    raw_level_out$COG_category <- category
    raw_level_out$taxon <- taxa_level
    all_output_raw[[combo_num]] <- raw_level_out
    combo_num <- combo_num + 1
  }
}

all_output <- do.call(rbind, all_output_raw)

all_output$fdr <- p.adjust(all_output$p, 'BH')

all_output$cog_descrip <- cog_descrip[all_output$category, 'V3']

write.table(x = all_output, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE,
            file = "/mfs/gdouglas/projects/ocean_hgt_zenodo/putative_hgt/cluster/clusterbased_posthoc_COG_enrich.tsv")
