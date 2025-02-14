rm(list = ls(all.names = TRUE))

hyperg_combined_info <- read.table("/mfs/gdouglas/projects/ocean_mags/networks/allsamples/clusterbased/metaG_hyperg_cooccur.allsamples.combined.tsv.gz",
                                   header=TRUE, sep="\t", stringsAsFactors = FALSE)

hyperg_combined_info$cooccur <- 'No'
hyperg_combined_info$cooccur_ratio <- hyperg_combined_info$cooccur_obs / hyperg_combined_info$cooccur_exp
hyperg_combined_info$cooccur[which(hyperg_combined_info$cooccur_BH < 0.05 & hyperg_combined_info$cooccur_ratio > 1)] <- 'Yes'

hyperg_combined_info$hgt <- 'No'
hyperg_combined_info$hgt[which(hyperg_combined_info$both_gene_count > 0)] <- 'Yes'

cor.test(hyperg_combined_info$cooccur_ratio, hyperg_combined_info$tip_dist, method = "spearman")

hyperg_combined_info_higheronly <- hyperg_combined_info[which(! hyperg_combined_info$diff_tax_level %in% c('Species', 'Strain')), ]

contingency_tab <- data.frame(matrix(NA, nrow=2, ncol=2))
colnames(contingency_tab) <- c('hgt_no', 'hgt_yes')
rownames(contingency_tab) <- c('cooccur_no', 'cooccur_yes')

contingency_tab['cooccur_no', 'hgt_no'] <- length(which(hyperg_combined_info_higheronly$hgt == 'No' & hyperg_combined_info_higheronly$cooccur == 'No'))
contingency_tab['cooccur_yes', 'hgt_no'] <- length(which(hyperg_combined_info_higheronly$hgt == 'No' & hyperg_combined_info_higheronly$cooccur == 'Yes'))
contingency_tab['cooccur_no', 'hgt_yes'] <- length(which(hyperg_combined_info_higheronly$hgt == 'Yes' & hyperg_combined_info_higheronly$cooccur == 'No'))
contingency_tab['cooccur_yes', 'hgt_yes'] <- length(which(hyperg_combined_info_higheronly$hgt == 'Yes' & hyperg_combined_info_higheronly$cooccur == 'Yes'))

fisher.test(contingency_tab)

cor.test(nocooccur_phylo_dist_higheronly$cooccur_ratio, nocooccur_phylo_dist_higheronly$tip_dist, method = "spearman")

hyperg_combined_info_loweronly <- hyperg_combined_info[which(hyperg_combined_info$diff_tax_level %in% c('Species', 'Strain')), ]
cor.test(hyperg_combined_info_loweronly$cooccur_ratio, hyperg_combined_info_loweronly$tip_dist, method = "spearman")

tipdist_hgt <- hyperg_combined_info_higheronly[which(hyperg_combined_info_higheronly$hgt == 'Yes'), 'tip_dist']
tipdist_nohgt <- hyperg_combined_info_higheronly[which(hyperg_combined_info_higheronly$hgt == 'No'), 'tip_dist']

mean(tipdist_hgt, na.rm = TRUE) - mean(tipdist_nohgt, na.rm = TRUE)
wilcox.test(tipdist_hgt, tipdist_nohgt, exact=FALSE)

tipdist_hgt_nocooccur <- hyperg_combined_info_higheronly[which(hyperg_combined_info_higheronly$hgt == 'Yes' & hyperg_combined_info_higheronly$cooccur == 'No'), 'tip_dist']
tipdist_hgt_cooccur <- hyperg_combined_info_higheronly[which(hyperg_combined_info_higheronly$hgt == 'Yes' & hyperg_combined_info_higheronly$cooccur == 'Yes'), 'tip_dist']

mean(tipdist_hgt_nocooccur, na.rm = TRUE) - mean(tipdist_hgt_cooccur, na.rm = TRUE)


