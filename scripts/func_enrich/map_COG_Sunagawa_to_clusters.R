rm(list = ls(all.names = TRUE))

# Get table of COGs + categories based on clusters.
COG2catgory_uniq <- readRDS("/mfs/gdouglas/db/COG_definitions/cog-20.to_category_collapse.rds")

gene_to_cluster <- read.table("/mfs/gdouglas/projects/ocean_mags/clusters/gene_to_cluster.tsv.gz",
                              header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = 1)

COG_annot <- read.table("/mfs/gdouglas/projects/ocean_mags/Sunagawa_dataset/genomes-cog-info.tsv.gz",
                        header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "", comment.char = "", row.names = 2)

# Drop unneeded columns.
COG_annot <- COG_annot[, c("all_COG", "COG_category")]

# Add in cluster as additional column.
COG_annot$cluster <- gene_to_cluster[rownames(COG_annot), "cluster_id"]

# Assess how much variation there is in COG category classification per cluster.
COG_annot_dedupped <- COG_annot
rownames(COG_annot_dedupped) <- NULL
COG_annot_dedupped <- COG_annot_dedupped[-which(duplicated(COG_annot_dedupped)), ]

nonconcordant_clusters <- unique(COG_annot_dedupped$cluster[which(duplicated(COG_annot_dedupped$cluster))])

# Use majority rule for discordant cases.
COG_annot_dedupped_concordant <- COG_annot_dedupped[which(! COG_annot_dedupped$cluster %in% nonconcordant_clusters), ]

COG_annot_discordant <- COG_annot[which(COG_annot$cluster %in% nonconcordant_clusters), ]

gene_to_cluster_discordant <- gene_to_cluster[which(gene_to_cluster$cluster_id %in% nonconcordant_clusters), , drop = FALSE]

gene_to_cluster_discordant_list <- collapse::rsplit(x = rownames(gene_to_cluster_discordant), fl = gene_to_cluster_discordant$cluster_id)

discordant_annot_raw <- parallel::mclapply(X = gene_to_cluster_discordant_list,

                                       FUN = function(members) {

                                         out_df <- data.frame(majority_COG = "", majority_COG_category = "", any_COG = "", any_COG_category = "")

                                         hits_annot <- COG_annot_discordant[which(rownames(COG_annot_discordant) %in% members), ]

                                         # Return NA only if no hits annotated.
                                         if (nrow(hits_annot) == 0) { return(out_df) }

                                         hit_annot_COG <- do.call(c, lapply(hits_annot$all_COG, function(x) { strsplit(x = x, split = ",")[[1]] }))

                                         # Annotate based on majority-rule if possible.
                                         majority_count_required <- ceiling((length(members) / 2) + 0.5)

                                         if (nrow(hits_annot) >= majority_count_required) {

                                           hit_annot_COG_tallies <- table(hit_annot_COG)

                                           majority_observed_i <- which(hit_annot_COG_tallies >= majority_count_required)

                                           if (length(majority_observed_i) > 0) {
                                             out_df$majority_COG <- paste(names(hit_annot_COG_tallies)[majority_observed_i], collapse = ",")
                                             raw_majority_COG_categories <- COG2catgory_uniq[names(hit_annot_COG_tallies)[majority_observed_i], "category"]
                                             out_df$majority_COG_category <- paste(unique(sort(do.call(c, lapply(raw_majority_COG_categories, function(x) { strsplit(x = x, split = ",")[[1]] })))), collapse = ",")
                                           }

                                         }


                                         # Also annotate based on "any" annotation approach.
                                         unique_COG_annot <- sort(unique(hit_annot_COG))
                                         out_df$any_COG <- paste(unique_COG_annot, collapse = ",")
                                         out_df$any_COG_category <- paste(unique(sort(do.call(c, lapply(COG2catgory_uniq[unique_COG_annot, "category"], function(x) { strsplit(x = x, split = ",")[[1]] })))), collapse = ",")

                                         return(out_df)
                                       },
                                       mc.cores = 100)


discordant_annot <- do.call(rbind, discordant_annot_raw)

discordant_annot <- discordant_annot[, c("majority_COG", "majority_COG_category")]
colnames(discordant_annot) <- c("all_COG", "COG_category")
discordant_annot$cluster <- rownames(discordant_annot)
rownames(discordant_annot) <- NULL

discordant_annot[which(discordant_annot$COG_category == ""), "COG_category"] <- "-"
discordant_annot[which(discordant_annot$all_COG == ""), "all_COG"] <- NA

cluster_annot <- rbind(COG_annot_dedupped_concordant, discordant_annot)

cluster_annot <- cluster_annot[, c("cluster", "all_COG", "COG_category")]

colnames(cluster_annot) <- c("cluster", "majority_rule_COG", "majority_rule_COG_category")

write.table(x = cluster_annot,
            file = "/mfs/gdouglas/projects/ocean_mags/Sunagawa_dataset/genomes-cdhit-cluster-cog-info.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
