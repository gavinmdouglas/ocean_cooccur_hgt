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

read_coverm_file <- function(file, subfolder, exp_rows = 15339) {

  cov_in <- read.table(file, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
  colnames(cov_in) <- c('genome', 'mean', 'rpkm', 'trimmed_mean', 'relabun', 'covered_bases', 'length', "readcount")

  if (nrow(cov_in) != exp_rows) { stop('CoverM file does not have expected number of rows!') }

  cov_in$breadth <- cov_in$covered_bases / cov_in$length
  orig_col <- colnames(cov_in)
  mgs_id <- gsub('\\..*$', '', basename(file))
  cov_in$mgs_sample <- mgs_id

  cov_in$subfolder <- subfolder
  cov_in <- cov_in[, c('subfolder', 'mgs_sample', orig_col)]

  return(cov_in)
}

coverm_read_by_folder <- function(folder) {
  all_files <- list.files(folder,
                          pattern = "cov.gz", full.names = TRUE)
  if (length(all_files) == 0) { stop("Error, no files")}
  filelist_raw <- mclapply(all_files, function(x) { read_coverm_file(file=x, subfolder=basename(folder)) }, mc.cores=64)
  return(data.frame(rbindlist(filelist_raw)))
}

write_gzip_table_w_rowcol <- function(x, outfile) {
  orig_col <- colnames(x)
  x$sample <- rownames(x)
  x <- x[, c("sample", orig_col)]

  gzfile_connection <-  gzfile(outfile, "w")

  write.table(x = x,
              file = gzfile_connection,
              sep = "\t",
              row.names = FALSE,
              col.names = TRUE,
              quote = FALSE)

  close(gzfile_connection)
}
