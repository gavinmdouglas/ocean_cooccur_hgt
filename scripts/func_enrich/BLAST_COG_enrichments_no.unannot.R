# Run COG category enrichments for 95-<99% and >=99% putatively transferred HGT hits at each taxonomic level.
# But this time, only consider COG-annotated clusters.

rm(list = ls(all.names = TRUE))

source('~/scripts/ocean_mag_hgt/scripts/func_enrich/enrich_function.R')

COG_category_descrip <- read.table("/mfs/gdouglas/db/COG_definitions/COG_category_descrip.tsv",
                                   header = FALSE, sep = "\t", row.names = 1)

representative_map <- read.table("/mfs/gdouglas/projects/ocean_mags/Sunagawa_dataset/gene-catalog-membership.tsv.gz",
                                 header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)

seq_lengths <- read.table("/mfs/gdouglas/projects/ocean_mags/Sunagawa_dataset/gene-catalog.fasta.lengths.tsv.gz",
                          header = TRUE, sep = "\t", row.names = 1)

cluster_annot <- read.table(file = "/mfs/gdouglas/projects/ocean_mags/Sunagawa_dataset/genomes-representative-cog-info.tsv.gz",
                            sep = "\t", row.names = 1, stringsAsFactors = FALSE, header = TRUE)

# Remove unannot clusters.
cluster_annot <- cluster_annot[which(cluster_annot$majority_rule_COG_category != "-"), ]

cluster_annot_mapping <- list()

background_all_genomes <- as.character()

for (COG_category in rownames(COG_category_descrip)) {
  cluster_annot_mapping[[COG_category]] <- rownames(cluster_annot)[grep(COG_category, cluster_annot$majority_rule_COG_category)]
  background_all_genomes <- unique(c(background_all_genomes, rownames(cluster_annot)[grep(COG_category, cluster_annot$majority_rule_COG_category)]))
}

to_ignore <- c("A", "B", "Y", "Z")

# Read in cross-taxa cluster info
clusters_crosstaxa_raw <- readRDS(file = "/mfs/gdouglas/projects/ocean_mags/blast_output/cluster_crosstaxa_summaries/2023.10.11_crosstaxa_clusters.rds")
names(clusters_crosstaxa_raw) <- gsub("pairwise_", "", names(clusters_crosstaxa_raw))

clusters_crosstaxa <- list()
for (category in names(clusters_crosstaxa_raw)) {
  clusters_crosstaxa[[category]] <- list()
  for (tax_level in names(clusters_crosstaxa_raw[[category]])) {
    clusters_crosstaxa[[category]][[tax_level]] <- clusters_crosstaxa_raw[[category]][[tax_level]]$clusters
    clusters_crosstaxa[[category]][[tax_level]] <- clusters_crosstaxa[[category]][[tax_level]][which(clusters_crosstaxa[[category]][[tax_level]] %in% rownames(cluster_annot))]
  }
}

all_output_raw <- list()

identities <- c('95', '99')
taxa_levels <- names(clusters_crosstaxa$`95`)

combo_num <- 1

for (identity in identities) {

  for (taxa_level in taxa_levels) {

    raw_level_out <- identify_enriched_categories(genes = clusters_crosstaxa[[identity]][[taxa_level]],
                                                  background = background_all_genomes[which(! background_all_genomes %in% clusters_crosstaxa[[identity]][[taxa_level]])],
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

write.table(x = all_output, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE,
            file = "/mfs/gdouglas/projects/ocean_mags/summary_files/COG_category_enrichment_no.unannot.tsv")
