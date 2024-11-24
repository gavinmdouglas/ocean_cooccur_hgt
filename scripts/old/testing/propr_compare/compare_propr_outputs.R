rm(list = ls(all.names = TRUE))

propr_compare <- read.table('~/tmp/propr_tab.tsv', header=FALSE, stringsAsFactors = FALSE, sep = " ")

cor.test(propr_compare$V2, propr_compare$V3, method="spearman")

plot(propr_compare$V2, propr_compare$V3)
