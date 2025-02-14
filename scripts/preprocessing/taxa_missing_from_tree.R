rm(list = ls(all.names = TRUE))

taxa <- read.table('/mfs/gdouglas/projects/ocean_mags/networks/allsamples/clusterbased/metaG_hyperg_cooccur.allsamples.combined.unique_taxa.txt',
                   stringsAsFactors = FALSE, header=FALSE)$V1

tree <- ape::read.tree('/mfs/gdouglas/projects/ocean_mags/phylogenetic_analyses/nico_tree/Aligned_SCGs.phy.treefile')

length(setdiff(taxa, tree$tip.label))

taxa_missing_from_tree <- setdiff(taxa, tree$tip.label)

write.table(x = taxa_missing_from_tree,
            file = '/mfs/gdouglas/projects/ocean_hgt_zenodo/taxa_subsets/taxa_missing_from_tree.txt',
            col.names = FALSE, row.names = FALSE, quote = FALSE)