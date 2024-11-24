rm(list = ls(all.names = TRUE))

hit_counts_per_region <- readRDS("/mfs/gdouglas/projects/ocean_mags/summary_files/per_hit_gene_counts.rds")

hit_counts_per_region <- hit_counts_per_region[which(! hit_counts_per_region$Level %in% c('Species', 'Strain')), ]

compute_length <- function(i) {

  hit_info <- strsplit(hit_counts_per_region$Hit[i], "\\|\\|")[[1]]

  start1 <- as.numeric(hit_info[2])
  stop1 <- as.numeric(hit_info[3])
  if (start1 < stop1) {
    length1 <- stop1 - start1 + 1
  } else {
    length1 <- start1 - stop1 + 1
  }

  start2 <- as.numeric(hit_info[5])
  stop2 <- as.numeric(hit_info[6])
  if (start2 < stop2) {
    length2 <- stop2 - start2 + 1
  } else {
    length2 <- start2 - stop2 + 1
  }

  return(mean(c(length1, length2)))
}

hit_lengths <- parallel::mclapply(1:nrow(hit_counts_per_region), compute_length, mc.cores=10)

hit_counts_per_region$length <- unlist(hit_lengths)

write.table(file = "/mfs/gdouglas/projects/ocean_hgt_zenodo/putative_hgt/blast/hit_gene_counts_and_lengths.tsv",
            x = hit_counts_per_region, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
