rm(list = ls(all.names = TRUE))

# Compute pairwise branch length distances between all tips in tree.

library(ape)
library(castor)

full_tree <- ape::read.tree('/mfs/nicot/Gavin/new_analysis/GToTree_New/Aligned_SCGs.phy.treefile')

tip_dist <- castor::get_all_pairwise_distances(tree = full_tree,
                                               only_clades = full_tree$tip.label)

tip_dist_df <- data.frame(as.matrix(tip_dist))

rownames(tip_dist_df) <- full_tree$tip.label
colnames(tip_dist_df) <- full_tree$tip.label

write.table(x = tip_dist_df,
            file = '/mfs/gdouglas/projects/ocean_mags/progenomes_analyses/genome_prep/high_qual_genome_tip_dist.tsv',
            col.names = NA,
            row.names = TRUE,
            quote = FALSE,
            sep = '\t')
