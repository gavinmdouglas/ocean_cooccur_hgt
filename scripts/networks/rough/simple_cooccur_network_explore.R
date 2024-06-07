rm(list = ls(all.names = TRUE))

library(data.table)
library(igraph)
library(reshape2)
library(parallel)
library(ggplot2)

dna_simple_cooccur_combined <- readRDS(file = "/mfs/gdouglas/projects/water_mags/coverm/network_working/dna_simple_cooccur.rds")

# Overall co-occurrence histogram (restricted to pairs with non-zero co-occurrence).
# Which is 
ggplot(data = dna_simple_cooccur_combined, aes(x = simple_cooccur)) +
  geom_histogram(bins=100) +
  theme_bw() +
  xlab("Overlap metric") +
  ylab("Frequency")

# Add tip distance and cross-taxon information to table.
tip_dist <- read.table("/mfs/gdouglas/projects/water_mags/water_mag_analysis/phylogenetic_analyses/tip_dist.tsv.gz",
                       header = TRUE, row.names = 1, sep = '\t', stringsAsFactors = FALSE)

taxa_breakdown <- read.table("/mfs/gdouglas/projects/water_mags/water_mag_analysis/mapfiles/MAG_taxa_breakdown.tsv.gz",
                             header = TRUE, row.names = 2, sep = '\t', stringsAsFactors = FALSE)

# Get basic centrality measures of graph.
dna_simple_cooccur_combined_graph <- igraph::graph_from_data_frame(d = dna_simple_cooccur_combined[, c("taxoni", "taxonj", "simple_cooccur")],
                                                                   directed = FALSE, vertices = NULL)

degree_breakdown <- igraph::degree(dna_simple_cooccur_combined_graph)
closeness_centrality <- igraph::closeness(dna_simple_cooccur_combined_graph, mode="all")
betweenness_centrality <- igraph::betweenness(dna_simple_cooccur_combined_graph, directed = FALSE)

degree_top10 <- names(sort(degree_breakdown, decreasing = TRUE))[1:10]
hist(degree_breakdown, main = "", xlab = "Number of degrees per node")
taxa_breakdown[degree_top10, ]

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
    
for (i in 1:1000) {
  taxoni <- as.character(dna_simple_cooccur_combined[i, "taxoni"])
  taxonj <- as.character(dna_simple_cooccur_combined[i, "taxonj"])
  tip_dist_by_pair <- c(tip_dist_by_pair, tip_dist[taxoni, taxonj])
  
  diff_tax_level <- "Identical"
  for (tax_level in tax_levels) {
    
    taxoni_tax <- taxa_breakdown[taxoni, tax_level]
    taxonj_tax <- taxa_breakdown[taxonj, tax_level]
  
    if (taxoni_tax != taxonj_tax) {
      diff_tax_level <- tax_level
      break
    }
  }
  tax_level_diff_by_pair <- c(tax_level_diff_by_pair, diff_tax_level)
}


# Co-occurrence vs. phylogenetic distance.

