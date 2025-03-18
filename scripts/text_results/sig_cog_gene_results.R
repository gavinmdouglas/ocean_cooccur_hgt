rm(list = ls(all.names = TRUE))

cog_results <- read.table('/mfs/gdouglas/projects/ocean_hgt_zenodo/putative_hgt/cluster/clusterbased_posthoc_COG_enrich.tsv.gz',
                          header=TRUE, sep ='\t', stringsAsFactors = FALSE, quote = '', comment.char = '')

cog_results[which(cog_results$fdr < 0.05), ]