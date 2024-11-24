rm(list = ls(all.names = TRUE))

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

