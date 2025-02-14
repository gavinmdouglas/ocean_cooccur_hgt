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

coverm_read_by_bioproject <- function(bioproject) {
  all_files <- list.files(paste0("/mfs/gdouglas/projects/ocean_mags/coverm/", bioproject),
                          pattern = "cov.gz", full.names = TRUE)
  if (length(all_files) == 0) { stop("Error, no files")}
  bioproject_raw <- mclapply(all_files, function(x) { read_coverm_file(file=x, bioproject=bioproject) }, mc.cores=64)

  return(data.frame(rbindlist(bioproject_raw)))
}

additional_coverm <- coverm_read_by_bioproject("additional_output")

PRJEB1787_coverm <- coverm_read_by_bioproject("PRJEB1787")
PRJEB1787_sample_run_map <- read.table("/mfs/gdouglas/projects/ocean_mags/coverm/combined_tables/Tara_dedup_id_mapping/PRJEB1787.tsv",
                                       header = TRUE, sep = "\t", row.names = 2, stringsAsFactors = FALSE)
PRJEB1787_coverm <- PRJEB1787_coverm[which(PRJEB1787_coverm$mgs_sample %in% rownames(PRJEB1787_sample_run_map)), ]
PRJEB1787_coverm$mgs_sample_recoded <- PRJEB1787_sample_run_map[PRJEB1787_coverm$mgs_sample, "sample_name"]
PRJEB1787_coverm$mgs_sample <- PRJEB1787_coverm$mgs_sample_recoded
PRJEB1787_coverm <- PRJEB1787_coverm[, -which(colnames(PRJEB1787_coverm) == "mgs_sample_recoded")]

PRJEB9740_coverm <- coverm_read_by_bioproject("PRJEB9740")
PRJEB9740_sample_run_map <- read.table("/mfs/gdouglas/projects/ocean_mags/coverm/combined_tables/Tara_dedup_id_mapping/PRJEB9740.tsv",
                                       header = TRUE, sep = "\t", row.names = 2, stringsAsFactors = FALSE)
PRJEB9740_coverm <- PRJEB9740_coverm[which(PRJEB9740_coverm$mgs_sample %in% rownames(PRJEB9740_sample_run_map)), ]
PRJEB9740_coverm$mgs_sample_recoded <- PRJEB9740_sample_run_map[PRJEB9740_coverm$mgs_sample, "sample_name"]
PRJEB9740_coverm$mgs_sample <- PRJEB9740_coverm$mgs_sample_recoded
PRJEB9740_coverm <- PRJEB9740_coverm[, -which(colnames(PRJEB9740_coverm) == "mgs_sample_recoded")]

combined_dna <- rbind(additional_coverm, PRJEB1787_coverm)
combined_dna <- rbind(combined_dna, PRJEB9740_coverm)
rownames(combined_dna) <- NULL

# Create presence tables.
dna_presence <- prep_presence_tab(cov_tab = combined_dna,
                                  presence_cutoff = 0.30)

# Create RPKM tables, based on the same samples/genomes as in the presence tables.
dna_rpkm <- reshape2::dcast(data = combined_dna,
                            formula = mgs_sample ~ genome,
                            value.var = "rpkm",
                            fill = 0)

dna_rpkm <- dna_rpkm[which(dna_rpkm$mgs_sample %in% rownames(dna_presence)),
                     c("mgs_sample", colnames(dna_presence))]

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

write_table_w_rowcol(x = dna_presence, outfile = "~/projects/ocean_mags/coverm/combined_tables/metaG_presence_allsamples.tsv")

write.table(x = dna_rpkm, file = "~/projects/ocean_mags/coverm/combined_tables/metaG_rpkm_allsamples.tsv",
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
