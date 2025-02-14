rm(list = ls(all.names = TRUE))

source('~/scripts/ocean_mag_hgt/scripts/func_enrich/enrich_function.R')

COG_category_descrip <- read.table("/mfs/gdouglas/db/COG_definitions/COG_category_descrip.tsv",
                                   header = FALSE, sep = "\t", row.names = 1)

func_info <- read.table('/mfs/gdouglas/projects/ocean_mags/functional_annot/gene_info_w_annot.tsv.gz',
                        stringsAsFactors = FALSE, sep = '\t', row.names=1, header=TRUE)

prepped_cooccur_hgt <- read.table('/mfs/gdouglas/projects/ocean_mags/networks/allsamples/clusterbased/metaG_hyperg_cooccur.allsamples.combined.hgt_only.tsv',
                                  sep = '\t', stringsAsFactors = FALSE, header=TRUE, row.names = 1)

cooccur_hgt_taxa_pairs <- rownames(prepped_cooccur_hgt)[which(prepped_cooccur_hgt$cooccur == 1)]
noncooccur_hgt_taxa_pairs <- rownames(prepped_cooccur_hgt)[which(prepped_cooccur_hgt$cooccur == 0)]

cluster_best_hit <- read.table('~/projects/ocean_hgt_zenodo/putative_hgt/cluster/all_best_hits.tsv.gz',
                               header=TRUE, sep='\t', stringsAsFactors = FALSE)

taxonomy <- read.table("/mfs/gdouglas/projects/ocean_hgt_zenodo/mapfiles/MAG_taxa_breakdown.tsv.gz",
                       header = TRUE, stringsAsFactors = FALSE, sep = "\t", row.names=1)

cluster_best_hit$taxon1 <- taxonomy[cluster_best_hit$gene1_genome, 'Taxon_ID']
cluster_best_hit$taxon2 <- taxonomy[cluster_best_hit$gene2_genome, 'Taxon_ID']

cluster_best_hit$taxa_pair <- NA
for (i in 1:nrow(cluster_best_hit)) {
 taxa_pair_raw <- as.character(cluster_best_hit[i, c('taxon1', 'taxon2')])
 cluster_best_hit$taxa_pair[i] <- paste(sort(taxa_pair_raw), collapse = ',')
}

cooccur_hgt_genes <- unique(c(cluster_best_hit[which(cluster_best_hit$taxa_pair %in% cooccur_hgt_taxa_pairs), 'gene1'],
                              cluster_best_hit[which(cluster_best_hit$taxa_pair %in% cooccur_hgt_taxa_pairs), 'gene2']))

noncooccur_hgt_genes <- unique(c(cluster_best_hit[which(cluster_best_hit$taxa_pair %in% noncooccur_hgt_taxa_pairs), 'gene1'],
                                 cluster_best_hit[which(cluster_best_hit$taxa_pair %in% noncooccur_hgt_taxa_pairs), 'gene2']))

func_info <- func_info[c(cooccur_hgt_genes, noncooccur_hgt_genes), ]

func_info_noCOGunclass <- func_info
func_info_noCOGunclass <- func_info_noCOGunclass[which(func_info_noCOGunclass$COG_category != '-'), ]

cluster_annot_mapping <- list()
for (COG_category in rownames(COG_category_descrip)) {
  cluster_annot_mapping[[COG_category]] <- rownames(func_info_noCOGunclass)[grep(COG_category, func_info_noCOGunclass$COG_category)]
}

COG_enrich <- identify_enriched_categories(genes = cooccur_hgt_genes,
                                           background = noncooccur_hgt_genes,
                                           gene_to_category_map = cluster_annot_mapping,
                                           to_ignore = c('A', 'B', 'Z'))
COG_categories_succinct <- read.table('/mfs/gdouglas/db/COG_definitions/COG_category_descrip_succinct.tsv',
                                      header = FALSE, sep = '\t', stringsAsFactors = FALSE, row.names = 1)

COG_enrich$COG_descrip <- COG_category_descrip[COG_enrich$category, 'V2']
COG_enrich$COG_descrip_succinct <- COG_categories_succinct[COG_enrich$category, 'V2']

write.table(x = COG_enrich,
            file = '/mfs/gdouglas/projects/ocean_hgt_zenodo/hgt_cooccur_enrich/cooccur_HGT_COG_enrich.tsv',
            sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)


# MGE-related annotation enrichment.
mge_map <- list()
mge_map[['proMGE']] <- rownames(func_info)[which(func_info$proMGE == 'Yes')]
mge_map[['scaffold_contains_provirus']] <- rownames(func_info)[which(func_info$scaffold_contains_provirus == 'Yes')]
mge_map[['Plasmid']] <- rownames(func_info)[which(func_info$genomad == 'Plasmid')]
mge_map[['Virus']] <- rownames(func_info)[which(func_info$genomad == 'Virus')]

mge_enrich <- identify_enriched_categories(genes = cooccur_hgt_genes,
                                           background = noncooccur_hgt_genes,
                                           gene_to_category_map = mge_map)

write.table(x = mge_enrich,
            file = '/mfs/gdouglas/projects/ocean_hgt_zenodo/hgt_cooccur_enrich/cooccur_HGT_MGE_enrich.tsv',
            sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)

