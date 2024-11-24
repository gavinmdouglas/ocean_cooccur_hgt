rm(list = ls(all.names = TRUE))

num_tested <- read.table('/mfs/gdouglas/projects/ocean_mags/species_DTL_analyses/num_tested_genes.txt',
                         header=TRUE, sep=' ', stringsAsFactors = FALSE)

dtl_summary <- read.table('/mfs/gdouglas/projects/ocean_mags/species_DTL_analyses/homer_rangerdtl_combined_summary.tsv',
                          header=TRUE, sep='\t', stringsAsFactors = FALSE, row.names = 1)

dtl_summary$num_tested <- NA

dtl_summary[num_tested$species, 'num_tested'] <- num_tested$num_tested_genes

prop_genes_w_transfer = dtl_summary$num_unique_genes_w_transfer / dtl_summary$num_tested

mean(prop_genes_w_transfer)
median(prop_genes_w_transfer)
sd(prop_genes_w_transfer)

sum(dtl_summary$num_total_gene_transfers)