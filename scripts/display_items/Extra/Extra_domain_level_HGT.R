# After running ~/scripts/ocean_mag_hgt/scripts/func_enrich/cluster.hgt_COG_enrichments_no.unannot.R (and not clearing memory).

COG_gene_descrip <- read.table("/mfs/gdouglas/db/COG_definitions/cog-20.def.tab",
                               header = FALSE, sep = "\t", row.names = 1, quote='', comment.char = '')

# Summary file for domain-level transfers.
domain_hit99_summary <- data.frame(matrix(NA, nrow = length(best_hits_by_tax[['99']]$Domain), ncol=5))
rownames(domain_hit99_summary) <- best_hits_by_tax[['99']]$Domain
colnames(domain_hit99_summary) <- c('cluster_id', 'identity_cutoff', 'COG_category', 'COG_gene', 'COG_gene_descrip')

domain_hit99_summary$cluster_id <- best_hits_by_tax[['99']]$Domain
domain_hit99_summary$identity_cutoff <- '>= 99%'

domain_hit95_summary <- data.frame(matrix(NA, nrow = length(best_hits_by_tax[['95']]$Domain), ncol=5))
rownames(domain_hit95_summary) <- best_hits_by_tax[['95']]$Domain
colnames(domain_hit95_summary) <- c('cluster_id', 'identity_cutoff', 'COG_category', 'COG_gene', 'COG_gene_descrip')

domain_hit95_summary$cluster_id <- best_hits_by_tax[['95']]$Domain
domain_hit95_summary$identity_cutoff <- '>= 95% and < 99%'

domain_hit_summary <- rbind(domain_hit99_summary, domain_hit95_summary)

hits_w_COG_cat <- intersect(rownames(domain_hit_summary), rownames(cluster_annot))

domain_hit_summary[hits_w_COG_cat, 'COG_category'] <- cluster_annot[hits_w_COG_cat, 'majority_rule_COG_category']
domain_hit_summary[hits_w_COG_cat, 'COG_gene'] <- cluster_annot[hits_w_COG_cat, 'majority_rule_COG']

for (i in 1:nrow(domain_hit_summary)) {

  if (is.na(domain_hit_summary[i, 'COG_gene'])) { next }

  cog_genes <- strsplit(x = domain_hit_summary[i, 'COG_gene'], split=',')[[1]]

  cog_descrip <- character()

  for (cog_gene in cog_genes) {
    if (! cog_gene %in% rownames(COG_gene_descrip)) {
      print("ERROR!")
      stop('ERROR')
    }
    cog_gene_info <- as.character(COG_gene_descrip[cog_gene, 2:6])
    cog_gene_info <- cog_gene_info[which(cog_gene_info != '')]
    cog_gene_info <- paste(cog_gene_info, collapse = ';')
    cog_descrip <- c(cog_descrip, cog_gene_info)
  }

  domain_hit_summary[i, 'COG_gene_descrip'] <- paste(cog_descrip, collapse = '|')
}

write.table(x = domain_hit_summary, file = '~/tmp/cross_domain_hits.tsv', sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)
