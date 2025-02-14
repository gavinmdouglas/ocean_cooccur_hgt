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

clusterbased_COG_enrich <- read.table('/mfs/gdouglas/projects/ocean_hgt_zenodo/hgt_cooccur_enrich/cooccur_HGT_COG_enrich.tsv.gz',
                                      header=TRUE, stringsAsFactors = FALSE, sep = '\t', row.names=1)


# Put all COG category odd's ratios in single table for plotting as heatmap.
COG_results <- data.frame(matrix(NA, nrow = 1, ncol = 20))
colnames(COG_results) <- c(expected_enriched, expected_depleted, expected_mixed)

COG_results[1, ] <- log2(clusterbased_COG_enrich[c(expected_enriched, expected_depleted, expected_mixed), 'OR'])
COG_results[which(clusterbased_COG_enrich[c(expected_enriched, expected_depleted, expected_mixed), 'BH'] >= 0.05)] <- NA

draw(Heatmap(matrix = as.matrix(COG_results),
                       heatmap_legend_param = list(title = expression("Sig. log"[2]*"(Odd's ratio)")),
                       na_col = 'grey95',
                       column_labels = COG_categories_succinct[colnames(COG_results), 'both'],
                       column_split = c(rep("Expected enriched", length(expected_enriched)),
                                        rep("Expected depleted", length(expected_depleted)),
                                        rep("Unclear expectation", length(expected_mixed))),
                       cluster_rows = FALSE,
                       cluster_columns = FALSE,

                       column_names_rot = 45,
                       column_gap = unit(5, "mm"),
                       col = colorRamp2(c(-4, -0.001, 0, 0.001, 4), c("slateblue1","lightsteelblue", "white", "lightpink", "red"))),
     padding = unit(c(2, 50, 2, 2), "mm"))

dev.off()
