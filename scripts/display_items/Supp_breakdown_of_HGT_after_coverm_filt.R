rm(list = ls(all.names = TRUE))

library(ComplexHeatmap)
library(gridGraphics)

# Basic table and visualization of % comparisons that were hits by taxonomic level and identity cut-off.
tax_levels <- c("Genus", "Family","Order", "Class", "Phylum", "Domain")

gene_info <- read.table("/mfs/gdouglas/projects/ocean_mags/Sunagawa_dataset/gene_info_allscaffolds.tsv.gz",
                        stringsAsFactors = FALSE, sep = '\t', header=TRUE, row.name=1)

gene_info_w_func_annot <- read.table("/mfs/gdouglas/projects/ocean_mags/functional_annot/gene_info_w_annot.tsv.gz",
                                     header=TRUE, sep="\t", stringsAsFactors = FALSE, row.names=1)

coverm_presence <- read.table('~/projects/ocean_mags/networks/combined_tables/metaG_presence_allsamples.tsv.gz',
                              header=TRUE, sep = '\t', stringsAsFactors = FALSE, row.names = 1)

taxonomy <- read.table('/mfs/gdouglas/projects/ocean_hgt_zenodo/mapfiles/MAG_taxa_breakdown.tsv.gz',
                       header=TRUE, sep = '\t', stringsAsFactors = FALSE, row.names=2)

present_genomes <- taxonomy[colnames(coverm_presence), 'MAG']

scaffolds_5000bp <- read.table('/mfs/gdouglas/projects/ocean_mags/Sunagawa_dataset/scaffolds_5000bp.txt',
                               stringsAsFactors = FALSE, header=FALSE)$V1

num_comparisons <- read.table('/mfs/gdouglas/projects/ocean_hgt_zenodo/putative_hgt/num_comparisons_per_inter.level.tsv.gz',
                              header = TRUE, sep = '\t', stringsAsFactors = FALSE, row.names = 1)

hgt_tab <- read.table('/mfs/gdouglas/projects/ocean_hgt_zenodo/putative_hgt/cluster/all_best_hits.tsv.gz',
                      header = TRUE, sep = '\t', stringsAsFactors = FALSE)

hgt_tab$gene1_scaffold <- gene_info[hgt_tab$gene1, 'scaffold']
hgt_tab$gene2_scaffold <- gene_info[hgt_tab$gene2, 'scaffold']

hgt_tab$scaffolds_retained <- FALSE

hgt_tab[which(hgt_tab$gene1_scaffold %in% scaffolds_5000bp & hgt_tab$gene2_scaffold %in% scaffolds_5000bp), 'scaffolds_retained'] <- TRUE

hgt_tab_genome_present <- hgt_tab
hgt_tab_genome_present <- hgt_tab_genome_present[-which(! hgt_tab_genome_present$gene1_genome %in% present_genomes), ]
hgt_tab_genome_present <- hgt_tab_genome_present[-which(! hgt_tab_genome_present$gene2_genome %in% present_genomes), ]

taxa_pairs <- c()
for (i in 1:nrow(hgt_tab_genome_present)) {
  taxa_pairs <- c(taxa_pairs, paste(sort(c(hgt_tab_genome_present[i, 'gene1_genome'], hgt_tab_genome_present[i, 'gene2_genome'])), collapse=','))
}
length(unique(taxa_pairs))
length(unique(c(hgt_tab_genome_present$gene1_genome, hgt_tab_genome_present$gene2_genome)))
# 715 pairwise genome combinations across 1,037 genomes based on genomes retained in coverm presence/absence table
# (i.e., those found at least 10 times).

taxa_pairs_orig <- c()
for (i in 1:nrow(hgt_tab)) {
  taxa_pairs_orig <- c(taxa_pairs_orig, paste(sort(c(hgt_tab[i, 'gene1_genome'], hgt_tab[i, 'gene2_genome'])), collapse=','))
}
length(unique(taxa_pairs_orig))
length(unique(c(hgt_tab$gene1_genome, hgt_tab$gene2_genome)))
# 2292 unique pairs and 2295 unique genomes.

# Now plot Phyla of genomes involved in HGT that were present before and after CoverM filter.
genome_summary <- read.table("/mfs/gdouglas/projects/ocean_mags/Sunagawa_dataset/genomes-summary.csv.gz",
                             sep = ",", stringsAsFactors = FALSE, header = TRUE, row.names = 1)

genome_summary$phylum <- genome_summary$GTDB.Taxonomy
genome_summary$phylum <- gsub(';c__.*$', '', genome_summary$phylum)

all_hgt_genomes <- unique(c(hgt_tab$gene1_genome, hgt_tab$gene2_genome))
hgt_genomes_coverm_retained <- unique(c(hgt_tab_genome_present$gene1_genome, hgt_tab_genome_present$gene2_genome))
hgt_genomes_coverm_filtered_out <- setdiff(all_hgt_genomes, hgt_genomes_coverm_retained)

map_tmp <-  read.table("/mfs/gdouglas/projects/ocean_hgt_zenodo/mapfiles/MAG_taxa_breakdown.tsv.gz",
                       header = TRUE, stringsAsFactors = FALSE, sep = "\t", row.names = 1)
analyzed_genomes <- rownames(map_tmp)

# Sanity check:
# identical(gsub(' ', '', map_tmp[analyzed_genomes, 'Phylum']), genome_summary[analyzed_genomes, 'phylum'])
# TRUE

unique_phyla <- unique(genome_summary$phylum)

phyla_breakdown <- data.frame(analyzed=rep(0, length(unique_phyla)),
                              hgt_all=rep(0, length(unique_phyla)),
                              hgt_coverm_retained=rep(0, length(unique_phyla)),
                              hgt_coverm_filtered_out=rep(0, length(unique_phyla)))
rownames(phyla_breakdown) <- unique_phyla

analyzed_summary <- genome_summary[analyzed_genomes, ]
analyzed_phyla_tally <- table(analyzed_summary$phylum)
for (phylum in names(analyzed_phyla_tally)) {
  phyla_breakdown[phylum, 'analyzed'] <- analyzed_phyla_tally[phylum]
}

hgt_all_summary <- genome_summary[all_hgt_genomes, ]
hgt_all_phyla_tally <- table(hgt_all_summary$phylum)
for (phylum in names(hgt_all_phyla_tally)) {
  phyla_breakdown[phylum, 'hgt_all'] <- hgt_all_phyla_tally[phylum]
}

hgt_coverm_retained_summary <- genome_summary[hgt_genomes_coverm_retained, ]
hgt_coverm_retained_phyla_tally <- table(hgt_coverm_retained_summary$phylum)
for (phylum in names(hgt_coverm_retained_phyla_tally)) {
  phyla_breakdown[phylum, 'hgt_coverm_retained'] <- hgt_coverm_retained_phyla_tally[phylum]
}

hgt_coverm_filtered_out_summary <- genome_summary[hgt_genomes_coverm_filtered_out, ]
hgt_coverm_filtered_out_phyla_tally <- table(hgt_coverm_filtered_out_summary$phylum)
for (phylum in names(hgt_coverm_filtered_out_phyla_tally)) {
  phyla_breakdown[phylum, 'hgt_coverm_filtered_out'] <- hgt_coverm_filtered_out_phyla_tally[phylum]
}

phyla_to_collapse_i <- which(rowSums(phyla_breakdown) < 15)
phyla_to_collapse_sums <- colSums(phyla_breakdown[phyla_to_collapse_i, ])
phyla_breakdown <- phyla_breakdown[-phyla_to_collapse_i, ]
phyla_breakdown <- phyla_breakdown[order(rowSums(phyla_breakdown), decreasing = TRUE), ]

phyla_breakdown['Rare (< 15)', ] <- phyla_to_collapse_sums

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

                                   column_labels = c('Analyzed', 'HGT all', 'HGT filt. out', 'HGT retained'),

                                   row_names_side = 'left',
                                   row_dend_side = 'right',
                                   column_names_rot = 45,

                                   cell_fun = function(j, i, x, y, width, height, fill) {
                                     if(! is.na(phyla_breakdown[i, j] > 0))
                                       grid.text(phyla_breakdown[i, j], x, y, gp = gpar(fontsize = 10), just = 'centre')
                                   })

phyla_breakdown_heatmap <- grid.grabExpr(draw(column_title = "", phyla_breakdown_heatmap))

ggsave(filename='/mfs/gdouglas/scripts/ocean_mag_hgt/display_items/Supp_phyla_breakdown_by_hgt_type.pdf',
       plot = phyla_breakdown_heatmap,
       height = 7, width = 6, dpi=600, device="pdf")
