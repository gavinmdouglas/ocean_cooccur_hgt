rm(list = ls(all.names = TRUE))

library(ComplexHeatmap)
library(gridGraphics)

genome_summary <- read.table("/mfs/gdouglas/projects/ocean_mags/Sunagawa_dataset/genomes-summary.csv.gz",
                             sep = ",", stringsAsFactors = FALSE, header = TRUE, row.names = 1)

genome_summary$phylum <- genome_summary$GTDB.Taxonomy
genome_summary$phylum <- gsub(';c__.*$', '', genome_summary$phylum)


analyzed_genomes <- read.table("/mfs/gdouglas/projects/ocean_hgt_zenodo/mapfiles/MAG_taxa_breakdown.tsv.gz",
                               header = TRUE, stringsAsFactors = FALSE, sep = "\t")$MAG

nonanalyzed_genomes <- setdiff(rownames(genome_summary), analyzed_genomes)



unique_phyla <- unique(genome_summary$phylum)

phyla_breakdown <- data.frame(analyzed=rep(0, length(unique_phyla)),
                              nonanalyzed=rep(0, length(unique_phyla)))
rownames(phyla_breakdown) <- unique_phyla

analyzed_summary <- genome_summary[analyzed_genomes, ]
analyzed_phyla_tally <- table(analyzed_summary$phylum)
for (phylum in names(analyzed_phyla_tally)) {
  phyla_breakdown[phylum, 'analyzed'] <- analyzed_phyla_tally[phylum]
}

nonanalyzed_summary <- genome_summary[nonanalyzed_genomes, ]
nonanalyzed_phyla_tally <- table(nonanalyzed_summary$phylum)
for (phylum in names(nonanalyzed_phyla_tally)) {
  phyla_breakdown[phylum, 'nonanalyzed'] <- nonanalyzed_phyla_tally[phylum]
}

phyla_to_collapse_i <- which(rowSums(phyla_breakdown) < 15)

phyla_to_collapse_sums <- colSums(phyla_breakdown[phyla_to_collapse_i, ])

phyla_NA <- phyla_breakdown['N/A', ]

phyla_breakdown <- phyla_breakdown[-phyla_to_collapse_i, ]
phyla_breakdown <- phyla_breakdown[-which(rownames(phyla_breakdown) == 'N/A'), ]

phyla_breakdown <- phyla_breakdown[order(rowSums(phyla_breakdown), decreasing = TRUE), ]


phyla_breakdown['Rare (< 15)', ] <- phyla_to_collapse_sums
phyla_breakdown['Unclassified', ] <- phyla_NA

bacteria_rows <- grep("d__Bacteria;", rownames(phyla_breakdown))
archaea_rows <- grep("d__Archaea;", rownames(phyla_breakdown))
other_rows <- setdiff(1:nrow(phyla_breakdown), c(bacteria_rows, archaea_rows))

rownames(phyla_breakdown) <- gsub('d__Bacteria;p__', '', rownames(phyla_breakdown))
rownames(phyla_breakdown) <- gsub('d__Archaea;p__', '', rownames(phyla_breakdown))

phyla_breakdown <- phyla_breakdown[c(bacteria_rows, archaea_rows, other_rows), ]

phyla_breakdown_percent <- as.matrix((phyla_breakdown / colSums(phyla_breakdown))) * 100
phyla_breakdown <- as.matrix(phyla_breakdown)

row_split_categories <- c(rep("Bacteria", length(bacteria_rows)),
                          rep("Archaea", length(archaea_rows)),
                          rep("Other", length(other_rows)))
row_split_categories <- factor(row_split_categories, levels=c("Bacteria", "Archaea", "Other"))

phyla_breakdown <- format(phyla_breakdown, big.mark = ',', justify = 'none', trim = TRUE)

phyla_breakdown_heatmap <- Heatmap(matrix = log10(phyla_breakdown_percent + 1),

                                    col = circlize::colorRamp2(c(0, 2), c('white', 'firebrick3')),

                                   heatmap_legend_param = list(title = expression(log[10](Percent + 1))),

                                    show_heatmap_legend = TRUE,

                                   row_split = row_split_categories,
                                   row_gap = unit(5, "mm"),
                                   row_title_rot=0,

                                    cluster_rows = FALSE,
                                    cluster_columns = FALSE,

                                   column_labels = c('Analyzed', 'Filtered out'),

                                    row_names_side = 'left',
                                    row_dend_side = 'right',
                                    column_names_rot = 45,

                                    cell_fun = function(j, i, x, y, width, height, fill) {
                                      if(! is.na(phyla_breakdown[i, j] > 0))
                                        grid.text(phyla_breakdown[i, j], x, y, gp = gpar(fontsize = 10), just = 'centre')
                                    })

phyla_breakdown_heatmap <- grid.grabExpr(draw(column_title = "", phyla_breakdown_heatmap))

ggsave(filename='/mfs/gdouglas/scripts/ocean_mag_hgt/display_items/Supp_phyla_breakdown.pdf',
       plot = phyla_breakdown_heatmap,
       height = 7, width = 5, dpi=600, device="pdf")
