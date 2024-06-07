rm(list = ls(all.names = TRUE))

library(reshape2)
library(parallel)
library(data.table)

prep_presence_tab <- function(cov_tab,
                              presence_cutoff,
                              min_samples_w_genome=10,
                              min_genomes_per_sample=1) {
  cov_tab$present <- 0
  
  cov_tab[which(cov_tab$breadth >= presence_cutoff), "present"] <- 1
  
  presence_matrix <- reshape2::dcast(data = cov_tab,
                                     formula = mgs_sample ~ genome,
                                     value.var = "present",
                                     fill = 0)
  
  rownames(presence_matrix) <- presence_matrix$mgs_sample
  presence_matrix <- presence_matrix[, which(colnames(presence_matrix) != "mgs_sample")]
  
  sample_tallies <- rowSums(presence_matrix)
  presence_matrix <- presence_matrix[which(sample_tallies >= min_genomes_per_sample), ]
  
  taxon_presence_tallies <- colSums(presence_matrix)
  presence_matrix <- presence_matrix[, which(taxon_presence_tallies >= min_samples_w_genome)]
  
  return(presence_matrix)
}

raw_in <- list()
read_coverm_file <- function(file, bioproject) {
  cov_in <- read.table(file, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
  
  if (bioproject == 'additional_output') {
    colnames(cov_in) <- c('genome', 'mean', 'rpkm', 'trimmed_mean', 'relabun', 'covered_bases', 'length', "readcount")
  } else {
    colnames(cov_in) <- c('genome', 'mean', 'rpkm', 'trimmed_mean', 'relabun', 'covered_bases', 'length')
    cov_in$readcount <- NA
  }
  
  cov_in$breadth <- cov_in$covered_bases / cov_in$length
  
  orig_col <- colnames(cov_in)
  
  mgs_id <- gsub('\\..*$', '', basename(file))
  cov_in$mgs_sample <- mgs_id
  
  cov_in$bioproject <- bioproject
  cov_in <- cov_in[, c('bioproject', 'mgs_sample', orig_col)]
  
  return(cov_in)
}

bioprojects <- c('additional_output')
bioproject_tallies <- data.frame(bioprojects = bioprojects,
                                 sample_number = NA)
rownames(bioproject_tallies) <- bioprojects

raw_in <- list()
for (bioproject in bioprojects) {
  all_files <- list.files(paste0("/mfs/gdouglas/projects/ocean_mags/coverm/", bioproject), 
                          pattern = "cov.gz", full.names = TRUE)
  bioproject_tallies[bioproject, 'sample_number'] <- length(all_files)
  bioproject_raw <- mclapply(all_files, function(x) { read_coverm_file(file=x, bioproject=bioproject) }, mc.cores=64)
  
  raw_in[[bioproject]] <- rbindlist(bioproject_raw)
}

combined_dna <- rbindlist(raw_in)
rownames(combined_dna) <- NULL
combined_dna$bioproject <- factor(combined_dna$bioproject, levels = bioprojects)

# Create presence tables.
dna_presence <- prep_presence_tab(cov_tab=combined_dna, presence_cutoff=0.30)

# Create RPKM tables, based on the same samples/genomes as in the presence tables.
dna_rpkm <- reshape2::dcast(data = combined_dna,
                            formula = mgs_sample ~ genome,
                            value.var = "rpkm",
                            fill = 0)

dna_rpkm <- dna_rpkm[which(dna_rpkm$mgs_sample %in% rownames(dna_presence)), c("mgs_sample", colnames(dna_presence))]

# Write output tables.
write_table_w_rowcol <- function(x, outfile) {
  orig_col <- colnames(x)
  x$sample <- rownames(x)
  x <- x[, c("sample", orig_col)]
  write.table(x = x,
              file = outfile,
              sep = "\t",
              row.names = FALSE,
              col.names = TRUE,
              quote = FALSE)
}

write_table_w_rowcol(x = dna_presence, outfile = "~/projects/ocean_mags/coverm/combined_tables/metaG_presence_extraonly.tsv")

write.table(x = dna_rpkm, file = "~/projects/ocean_mags/coverm/combined_tables/metaG_rpkm_extraonly.tsv",
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
