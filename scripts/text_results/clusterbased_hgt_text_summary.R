rm(list = ls(all.names = TRUE))

clusterbased_hits <- read.table('/mfs/gdouglas/projects/ocean_hgt_zenodo/putative_hgt/cluster/cross_level_tallies_norm.tsv.gz',
                                stringsAsFactors = FALSE, sep = '\t', header = TRUE)

hgt_tab <- read.table('/mfs/gdouglas/projects/ocean_hgt_zenodo/putative_hgt/cluster/all_best_hits.tsv.gz',
                      header = TRUE, sep = '\t', stringsAsFactors = FALSE)

# Total number
sum(clusterbased_hits$both)
nrow(hgt_tab)


# Number of HGT events per genome pair (of those with at least one)
taxa_pairs <- c()
for (i in 1:nrow(hgt_tab)) {
  taxa_pairs <- c(taxa_pairs, paste(sort(c(hgt_tab[i, 'gene1_genome'], hgt_tab[i, 'gene2_genome'])), collapse=','))
}

median((table(taxa_pairs)))
mean(table(taxa_pairs))



# And sanity check on Supplementary table 1:
map_tmp <-  read.table("/mfs/gdouglas/projects/ocean_hgt_zenodo/mapfiles/MAG_taxa_breakdown.tsv.gz",
                       header = TRUE, stringsAsFactors = FALSE, sep = "\t", row.names = 1)


table(map_tmp[unique(c(hgt_tab$gene1_genome, hgt_tab$gene2_genome)), 'Phylum'])

coverm_presence <- read.table('~/projects/ocean_mags/networks/combined_tables/metaG_presence_allsamples.tsv.gz',
                              header=TRUE, sep = '\t', stringsAsFactors = FALSE, row.names = 1)


hgt_tab$taxon1 <- map_tmp[hgt_tab$gene1_genome, 'Taxon_ID']
hgt_tab$taxon2 <- map_tmp[hgt_tab$gene2_genome, 'Taxon_ID']

hgt_tab_filt <- hgt_tab[hgt_tab$taxon1 %in% colnames(coverm_presence) & hgt_tab$taxon2 %in% colnames(coverm_presence), ]


table(map_tmp[unique(c(hgt_tab_filt$gene1_genome, hgt_tab_filt$gene2_genome)), 'Phylum'])

# Percent of genome excluded at CoverM stage.

(1 - (ncol(coverm_presence)/15339)) * 100

# Percent of HGT genomes lost.
num_orig_hgt_genomes <- length(unique(c(hgt_tab$gene1_genome, hgt_tab$gene2_genome)))
num_filt_hgt_genomes <- length(unique(c(hgt_tab_filt$gene1_genome, hgt_tab_filt$gene2_genome)))

(1 - (num_filt_hgt_genomes / num_orig_hgt_genomes)) * 100
