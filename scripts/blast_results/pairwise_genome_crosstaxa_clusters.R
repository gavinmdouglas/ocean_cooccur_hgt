rm(list = ls(all.names = TRUE))

gene_info <- read.table('/data2/gdouglas/projects/water_mag_analysis/Sunagawa_dataset/gene-catalog-membership.tsv.gz',
                        header = TRUE, sep = '\t', row.names = 1, stringsAsFactors = FALSE)

identities <- c('95', '99')

cluster_crosstaxa_pairwise_approach <- list()

for (identity in identities) {
  
  category_name <- paste('pairwise', identity, sep = '_')
  
  cluster_crosstaxa_pairwise_approach[[category_name]] <- list()
  
  for (tax_level in c('Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'Strain')) {
    
    bed_path <- paste('/data2/gdouglas/projects/water_mag_analysis/blast_output/intersections/blast_hit_gene_beds/Hits_',
                      tax_level,
                      '_',
                      identity,
                      '_genes.bed',
                      sep = '')
    
    wc_cmd <- paste('wc -l', bed_path)
    
    if (as.integer(strsplit(x = as.character(system(wc_cmd, intern = TRUE)), split = ' ')[[1]][1]) != 0) {
      
      gene_hits <- read.table(bed_path, sep = '\t', stringsAsFactors = FALSE, quote = '', comment.char = '')$V4
      
      gene_hits <- unique(gene_hits[which(gene_hits %in% rownames(gene_info))])
      
      cluster_hits <- unique(gene_info[gene_hits, 'representative'])
      
      cluster_crosstaxa_pairwise_approach[[category_name]][[tax_level]] <- list('clusters' = cluster_hits,
                                                                                'genes' = gene_hits)
      
    } else {
      
      cluster_crosstaxa_pairwise_approach[[category_name]][[tax_level]] <- list('clusters' = character(),
                                                                                'genes' = character())
      
    }
    
  }
  
}

RDS_filename <- '/data2/gdouglas/projects/water_mag_analysis/blast_output/cluster_crosstaxa_summaries/2023.10.11_crosstaxa_clusters.rds'

saveRDS(object = cluster_crosstaxa_pairwise_approach,
        file = RDS_filename)
