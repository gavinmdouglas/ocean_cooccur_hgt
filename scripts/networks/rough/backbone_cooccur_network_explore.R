rm(list = ls(all.names = TRUE))

library(data.table)
library(igraph)
library(reshape2)
library(parallel)
library(ggplot2)

dna_backbone_cooccur <- readRDS(file = "/mfs/gdouglas/projects/water_mags/coverm/network_working/dna_backbone_cooccur.rds")

# Get basic centrality measures of graph.
dna_backbone_cooccur_graph <- igraph::graph_from_data_frame(d = dna_backbone_cooccur[, c("taxon_i", "taxon_j", "P")],
                                                            directed = FALSE, vertices = NULL)

plot(dna_backbone_cooccur_graph, vertex.label=NA, vertex.size=0.5)

taxa_breakdown <- read.table("/mfs/gdouglas/projects/water_mags/water_mag_analysis/mapfiles/MAG_taxa_breakdown.tsv.gz",
                             header = TRUE, row.names = 2, sep = '\t', stringsAsFactors = FALSE)

degree_breakdown <- igraph::degree(dna_backbone_cooccur_graph)
degree_top10 <- names(sort(degree_breakdown, decreasing = TRUE))[1:10]
taxa_breakdown[degree_top10, ]

hist(degree_breakdown, main = "", xlab = "Number of degrees per node")


closeness_centrality <- igraph::closeness(dna_backbone_cooccur_graph, mode="all")
betweenness_centrality <- igraph::betweenness(dna_backbone_cooccur_graph, directed = FALSE)


dna_backbone_cooccur_sorted <- dna_backbone_cooccur[order(dna_backbone_cooccur$P, decreasing=FALSE), ]
dna_backbone_cooccur_top5000 <- dna_backbone_cooccur_sorted[1:5000, ]

dna_backbone_cooccur_top5000_graph <- igraph::graph_from_data_frame(d = dna_backbone_cooccur_top5000[, c("taxon_i", "taxon_j", "P")],
                                                                    directed = FALSE, vertices = NULL)

plot(dna_backbone_cooccur_top5000_graph, vertex.label=NA, vertex.size=0.5)

taxa_breakdown <- read.table("/mfs/gdouglas/projects/water_mags/water_mag_analysis/mapfiles/MAG_taxa_breakdown.tsv.gz",
                             header = TRUE, row.names = 2, sep = '\t', stringsAsFactors = FALSE)

degree_breakdown <- igraph::degree(dna_backbone_cooccur_graph)
degree_top10 <- names(sort(degree_breakdown, decreasing = TRUE))[1:10]
taxa_breakdown[degree_top10, ]

hist(degree_breakdown, main = "", xlab = "Number of degrees per node")


closeness_centrality <- igraph::closeness(dna_backbone_cooccur_graph, mode="all")
betweenness_centrality <- igraph::betweenness(dna_backbone_cooccur_graph, directed = FALSE)


# Add tip distance and cross-taxon information to table.
tip_dist <- read.table("/mfs/gdouglas/projects/water_mags/water_mag_analysis/phylogenetic_analyses/tip_dist.tsv.gz",
                       header = TRUE, row.names = 1, sep = '\t', stringsAsFactors = FALSE)


tip_dist_by_pair <- numeric()

tax_levels <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain")
tax_level_diff_by_pair <- character()

tax_info_raw <- mclapply(1:nrow(dna_simple_cooccur_combined),
                         function(i) {
                           taxoni <- as.character(dna_simple_cooccur_combined[i, "taxoni"])
                           taxonj <- as.character(dna_simple_cooccur_combined[i, "taxonj"])
                           tip_dist_by_pair <- tip_dist[taxoni, taxonj]
                           
                           diff_tax_level <- "Identical"
                           for (tax_level in tax_levels) {
                             
                             taxoni_tax <- taxa_breakdown[taxoni, tax_level]
                             taxonj_tax <- taxa_breakdown[taxonj, tax_level]
                             
                             if (taxoni_tax != taxonj_tax) {
                               diff_tax_level <- tax_level
                               break
                             }
                           }
                           
                           return(list(tip_dist_by_pair=tip_dist_by_pair, diff_tax_level=diff_tax_level))
                           
                         },
                         mc.cores=100)
