rm(list = ls(all.names = TRUE))

cluster_best_hit <- read.table('~/projects/ocean_hgt_zenodo/putative_hgt/cluster/all_best_hits.tsv.gz',
                               header=TRUE, sep='\t', stringsAsFactors = FALSE)

func_info <- read.table('/mfs/gdouglas/projects/ocean_mags/functional_annot/gene_info_w_annot.tsv.gz',
                        stringsAsFactors = FALSE, sep = '\t', row.names=1, header=TRUE)

hgt_genes <- unique(c(cluster_best_hit$gene1, cluster_best_hit$gene2))
nonhgt_genes <- setdiff(rownames(func_info), hgt_genes)

hgt_genes_proMGE <- table(func_info[hgt_genes, 'proMGE'])
nonhgt_genes_proMGE <- table(func_info[nonhgt_genes, 'proMGE'])
(hgt_genes_proMGE['Yes'] / sum(hgt_genes_proMGE)) / (nonhgt_genes_proMGE['Yes'] / sum(nonhgt_genes_proMGE))

hgt_genes_scaffold_contains_provirus <- table(func_info[hgt_genes, 'scaffold_contains_provirus'])
nonhgt_genes_scaffold_contains_provirus <- table(func_info[nonhgt_genes, 'scaffold_contains_provirus'])
(hgt_genes_scaffold_contains_provirus['Yes'] / sum(hgt_genes_scaffold_contains_provirus)) / (nonhgt_genes_scaffold_contains_provirus['Yes'] / sum(nonhgt_genes_scaffold_contains_provirus))

hgt_genes_genomad <- table(func_info[hgt_genes, 'genomad'])
nonhgt_genes_genomad <- table(func_info[nonhgt_genes, 'genomad'])

(hgt_genes_genomad['Plasmid'] / sum(hgt_genes_genomad)) / (nonhgt_genes_genomad['Plasmid'] / sum(nonhgt_genes_genomad))
(hgt_genes_genomad['Virus'] / sum(hgt_genes_genomad)) / (nonhgt_genes_genomad['Virus'] / sum(nonhgt_genes_genomad))

