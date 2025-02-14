rm(list = ls(all.names = TRUE))

library(circlize)
library(ComplexHeatmap)
library(gridtext)

COG_categories_succinct <- read.table('/mfs/gdouglas/db/COG_definitions/COG_category_descrip_succinct.tsv',
                                      header = FALSE, sep = '\t', stringsAsFactors = FALSE, row.names = 1)

COG_categories_succinct$both <- paste(COG_categories_succinct$V2, ' (', rownames(COG_categories_succinct), ')', sep = '')

expected_hgt <- read.table('/mfs/gdouglas/projects/ocean_hgt_zenodo/mapfiles/Dmitrijeva2024_COG_category_HGT_summary.tsv.gz',
                           header = TRUE, sep = '\t', stringsAsFactors = FALSE)

expected_enriched <- expected_hgt$COG_category[which(expected_hgt$Custom_summary == 'Enriched')]
expected_enriched <- expected_enriched[which(expected_enriched != "K (reg.)")]

expected_depleted <- expected_hgt$COG_category[which(expected_hgt$Custom_summary == 'Depleted')]

expected_mixed <- expected_hgt$COG_category[which(expected_hgt$Custom_summary == 'Mixed')]
expected_mixed[expected_mixed == "K (all)"] <- "K"

clusterbased_COG_enrich <- read.table('/mfs/gdouglas/projects/ocean_hgt_zenodo/putative_hgt/cluster/clusterbased_COG_category_enrichment_no.unannot.tsv.gz',
                                      header=TRUE, stringsAsFactors = FALSE, sep = '\t')


# Put all COG category odd's ratios in single table for plotting as heatmap.
COG_results <- data.frame(matrix(NA, nrow = 12, ncol = 22))
colnames(COG_results) <- c("Approach", "Taxon_level", expected_enriched, expected_depleted, expected_mixed)

COG_results$Approach <- c(rep("HGT\n>= 95%\nand\n< 99%", 6), rep("HGT\n>= 99%", 6))
COG_results$Taxon_level <- rep(c("Genus", "Family", "Order", "Class", "Phylum", "Domain"), 2)

# Also separately take a look at the cluster-based "both" results to make sure they make sense, but no need to plot as well.
cluster_COG_both_results <- data.frame(matrix(NA, nrow = 6, ncol = 21))
colnames(cluster_COG_both_results) <- c("Taxon_level", expected_enriched, expected_depleted, expected_mixed)
cluster_COG_both_results$Taxon_level <- c("Genus", "Family", "Order", "Class", "Phylum", "Domain")

for (i in 1:nrow(clusterbased_COG_enrich)) {
 identity_cutoff <- clusterbased_COG_enrich[i, "identity"]
 taxon_level <- clusterbased_COG_enrich[i, "taxon"]
 COG_category <- clusterbased_COG_enrich[i, "category"]

 if (! COG_category %in% colnames(COG_results)) { next }

 if (clusterbased_COG_enrich[i, "fdr"] >= 0.05) { next }

 if (identity_cutoff == "95") {
   COG_results[which(COG_results$Approach == "HGT\n>= 95%\nand\n< 99%" & COG_results$Taxon_level == taxon_level), COG_category] <- log2(clusterbased_COG_enrich[i, "OR"])
 } else if (identity_cutoff == "99") {
   COG_results[which(COG_results$Approach == "HGT\n>= 99%" & COG_results$Taxon_level == taxon_level), COG_category] <- log2(clusterbased_COG_enrich[i, "OR"])
 } else if (identity_cutoff == "both") {
   cluster_COG_both_results[which(cluster_COG_both_results$Taxon_level == taxon_level), COG_category] <- log2(clusterbased_COG_enrich[i, "OR"])
 }
}

COG_results_ORs <- as.matrix(COG_results[, 3:ncol(COG_results)])


COG_results$Approach <- gsub(" \\(", "\n\\(", COG_results$Approach)
COG_results$Approach <- factor(COG_results$Approach,
                               levels = c("HGT\n>= 95%\nand\n< 99%",
                                          "HGT\n>= 99%"))


pdf("/mfs/gdouglas/scripts/ocean_mag_hgt/display_items/Main_COG_enrich_heatmap.pdf", width = 10, height = 6.5)

draw(Heatmap(matrix = COG_results_ORs,
                       heatmap_legend_param = list(title = expression("Sig. log"[2]*"(Odd's ratio)")),
                       row_names_side = "left",
                       row_labels = COG_results$Taxon_level,
                       na_col = 'grey95',
                       column_labels = COG_categories_succinct[colnames(COG_results_ORs), 'both'],
                       column_split = c(rep("Expected enriched", length(expected_enriched)),
                                        rep("Expected depleted", length(expected_depleted)),
                                        rep("Unclear expectation", length(expected_mixed))),
                       row_split = COG_results$Approach,
                       row_gap = unit(5, "mm"),
                       cluster_rows = FALSE,
                       cluster_columns = FALSE,
                       row_title_rot = 0,

                       column_names_rot = 45,
                       column_gap = unit(5, "mm"),
                       col = colorRamp2(c(-4, -0.001, 0, 0.001, 4), c("slateblue1","lightsteelblue", "white", "lightpink", "red"))),
     padding = unit(c(2, 20, 2, 2), "mm"))

dev.off()
