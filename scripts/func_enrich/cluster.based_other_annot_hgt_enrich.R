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

proMGE_tally_tab <- matrix(c(as.integer(nonhgt_genes_proMGE), as.integer(hgt_genes_proMGE)), nrow=2, ncol=2)
fisher.test(proMGE_tally_tab)


hgt_genes_scaffold_contains_provirus <- table(func_info[hgt_genes, 'scaffold_contains_provirus'])
nonhgt_genes_scaffold_contains_provirus <- table(func_info[nonhgt_genes, 'scaffold_contains_provirus'])
(hgt_genes_scaffold_contains_provirus['Yes'] / sum(hgt_genes_scaffold_contains_provirus)) / (nonhgt_genes_scaffold_contains_provirus['Yes'] / sum(nonhgt_genes_scaffold_contains_provirus))

scaffold_contains_provirus_tally_tab <- matrix(c(as.integer(nonhgt_genes_scaffold_contains_provirus), as.integer(hgt_genes_scaffold_contains_provirus)), nrow=2, ncol=2)
fisher.test(scaffold_contains_provirus_tally_tab)


hgt_genes_genomad <- table(func_info[hgt_genes, 'genomad'])
nonhgt_genes_genomad <- table(func_info[nonhgt_genes, 'genomad'])

(hgt_genes_genomad['Plasmid'] / sum(hgt_genes_genomad)) / (nonhgt_genes_genomad['Plasmid'] / sum(nonhgt_genes_genomad))
(hgt_genes_genomad['Virus'] / sum(hgt_genes_genomad)) / (nonhgt_genes_genomad['Virus'] / sum(nonhgt_genes_genomad))

genomad_virus_tally_tab <- matrix(c(nonhgt_genes_genomad['Other'] + nonhgt_genes_genomad['Plasmid'], hgt_genes_genomad['Other'] + hgt_genes_genomad['Plasmid'],
                                    nonhgt_genes_genomad['Virus'], hgt_genes_genomad['Virus']), nrow=2, ncol=2)
fisher.test(genomad_virus_tally_tab)


genomad_plasmid_tally_tab <- matrix(c(nonhgt_genes_genomad['Other'] + nonhgt_genes_genomad['Virus'], hgt_genes_genomad['Other'] + hgt_genes_genomad['Virus'],
                                    nonhgt_genes_genomad['Plasmid'], hgt_genes_genomad['Plasmid']), nrow=2, ncol=2)
fisher.test(genomad_plasmid_tally_tab)