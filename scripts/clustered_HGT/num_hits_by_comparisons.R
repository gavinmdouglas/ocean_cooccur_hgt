rm(list = ls(all.names = TRUE))

# Basic table and visualization of % comparisons that were hits by taxonomic level and identity cut-off.
tax_levels <- c("Genus", "Family","Order", "Class", "Phylum", "Domain")

num_comparisons <- read.table('/mfs/gdouglas/projects/ocean_hgt_zenodo/putative_hgt/num_comparisons_per_inter.level.tsv.gz',
                              header = TRUE, sep = '\t', stringsAsFactors = FALSE, row.names = 1)

hgt_tab <- read.table('/mfs/gdouglas/projects/ocean_hgt_zenodo/putative_hgt/cluster/all_best_hits.tsv.gz',
                      header = TRUE, sep = '\t', stringsAsFactors = FALSE)

hist(hgt_tab$identity, main="Cluster-based HGT", breaks=50, xlab="Percent identity")

hgt_tab_95 <- hgt_tab[which(hgt_tab$identity >= 95.0 & hgt_tab$identity < 99.0), ]
hgt_tab_99 <- hgt_tab[which(hgt_tab$identity >= 99.0), ]

num_hits <- data.frame(Level=tax_levels,
                        num_comparisons=NA,
                        cluster_95_99=NA,
                        cluster_99=NA)
rownames(num_hits) <- num_hits$Level
num_hits$num_comparisons <- num_comparisons[rownames(num_hits), 1]

for (tax_level in tax_levels) {
  num_hits[tax_level, "cluster_95_99"] <- length(which(hgt_tab_95$highest_tax_diff == tax_level))
  num_hits[tax_level, "cluster_99"] <- length(which(hgt_tab_99$highest_tax_diff == tax_level))
}

num_hits$both <- num_hits$cluster_95_99 + num_hits$cluster_99

num_hits$cluster_95_99.norm <- num_hits$cluster_95_99 / num_hits$num_comparisons
num_hits$cluster_99.norm <- num_hits$cluster_99 / num_hits$num_comparisons
num_hits$both.norm <- num_hits$both / num_hits$num_comparisons

write.table(x = num_hits,
            file = "/mfs/gdouglas/projects/ocean_hgt_zenodo/putative_hgt/cluster/cross_level_tallies_norm.tsv",
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
