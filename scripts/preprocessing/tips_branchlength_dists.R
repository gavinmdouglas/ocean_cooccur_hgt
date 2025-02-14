rm(list = ls(all.names = TRUE))

# Compute pairwise branch length distances between all tips in tree.

library(ape)
library(castor)

full_tree <- ape::read.tree('/data2/gdouglas/projects/water_mag_analysis/phylogenetic_analyses/nico_tree/Aligned_SCGs.phy.treefile')

tip_dist <- castor::get_all_pairwise_distances(tree = full_tree,
                                               only_clades = full_tree$tip.label)

tip_dist_df <- data.frame(as.matrix(tip_dist))

rownames(tip_dist_df) <- full_tree$tip.label
colnames(tip_dist_df) <- full_tree$tip.label

write.table(x = tip_dist_df,
            file = '/data2/gdouglas/projects/water_mag_analysis/phylogenetic_analyses/tip_dist.tsv',
            col.names = NA,
            row.names = TRUE,
            quote = FALSE,
            sep = '\t')
