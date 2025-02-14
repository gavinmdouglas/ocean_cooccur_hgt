rm(list = ls(all.names = TRUE))

taxonomy_raw <- readRDS("/mfs/gdouglas/projects/ocean_mags/progenomes_analyses/genome_prep/taxid_taxize_download.rds")

taxids <- read.table('/mfs/gdouglas/projects/ocean_mags/progenomes_analyses/genome_prep/all_taxids.txt', stringsAsFactors = FALSE, header=FALSE)$V1

taxa_levels <- c('superkingdom', 'kingdom', 'phylum', 'subphylum', 'superclass', 'class', 'subclass',
                 'infraclass', 'superorder', 'order', 'suborder', 'infraorder', 'parvorder', 'superfamily',
                 'family', 'subfamily', 'genus', 'species', 'subspecies', 'strain')

taxa <- data.frame(matrix(NA, nrow = length(taxids), ncol = length(taxa_levels)))
rownames(taxa) <- taxids
colnames(taxa) <- taxa_levels

taxonomy <- list()
all_ranks <- character()
for (i in 1:length(taxonomy_raw)) {
  i_name <- names(taxonomy_raw[[i]])
  taxonomy[[i_name]] <- data.frame(taxonomy_raw[[i]][[i_name]])
  taxonomy[[i_name]] <- taxonomy[[i_name]][which(taxonomy[[i_name]]$rank != 'clade'), ]
  taxonomy[[i_name]] <- taxonomy[[i_name]][which(taxonomy[[i_name]]$rank != 'no rank'), ]
  all_ranks <- c(all_ranks, taxonomy[[i_name]]$rank)
  for (taxa_level in taxa_levels) {
    if (taxa_level %in% taxonomy[[i_name]]$rank) {
      rank_name <- taxonomy[[i_name]][which(taxonomy[[i_name]]$rank == taxa_level), 'name']
      if (length(rank_name) != 1) {
        stop('ERROR')
      } else {
        taxa[i_name, taxa_level] <- rank_name
      }
    }
  }
}

clean_tab <- data.frame(taxid = rownames(taxa),
                        Domain = taxa$superkingdom,
                        Phylum = paste(taxa$superkingdom, taxa$phylum, sep = '; '),
                        Class = paste(taxa$superkingdom, taxa$phylum, taxa$class, sep = '; '),
                        Order = paste(taxa$superkingdom, taxa$phylum, taxa$class, taxa$order, sep = '; '),
                        Family = paste(taxa$superkingdom, taxa$phylum, taxa$class, taxa$order, taxa$family, sep = '; '),
                        Genus = paste(taxa$superkingdom, taxa$phylum, taxa$class, taxa$order, taxa$family, taxa$genus, sep = '; '),
                        Species = paste(taxa$superkingdom, taxa$phylum, taxa$class, taxa$order, taxa$family, taxa$genus, taxa$species, sep = '; '),
                        Strain = paste(taxa$superkingdom, taxa$phylum, taxa$class, taxa$order, taxa$family, taxa$genus, taxa$species, taxa$strain, sep = '; '))

write.table(file = gzfile('/mfs/gdouglas/projects/ocean_mags/progenomes_analyses/genome_prep/progenomes_ncbi_taxonomy.tsv.gz'),
            x = clean_tab, col.names = TRUE, row.names = FALSE, sep = '\t', quote=FALSE)
