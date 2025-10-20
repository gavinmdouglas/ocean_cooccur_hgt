rm(list = ls(all.names = TRUE))

hit_counts_per_region <- read.table("/mfs/gdouglas/projects/ocean_hgt_zenodo/putative_hgt/blast/hit_gene_counts_and_lengths.tsv.gz",
                                   header=TRUE, sep = "\t", stringsAsFactors = FALSE)

mean(hit_counts_per_region[hit_counts_per_region$Identity == "Identity >= 99%", "Num_covered_genes"])
sd(hit_counts_per_region[hit_counts_per_region$Identity == "Identity >= 99%", "Num_covered_genes"])

mean(hit_counts_per_region[hit_counts_per_region$Identity == "Identity >= 95% and < 99%", "Num_covered_genes"])
sd(hit_counts_per_region[hit_counts_per_region$Identity == "Identity >= 95% and < 99%", "Num_covered_genes"])
