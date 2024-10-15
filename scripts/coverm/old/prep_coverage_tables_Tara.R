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

bioprojects <- c('PRJEB1787', 'PRJEB9740', 'PRJEB6608')

bioproject_tallies <- data.frame(bioprojects = bioprojects,
                                 sample_number = NA)
rownames(bioproject_tallies) <- bioprojects

raw_in <- list()

for (bioproject in bioprojects) {
  
  all_files <- list.files(paste0("~/projects/ocean_mags/coverm/", bioproject), 
                          pattern = "cov.gz", full.names = TRUE)
  
  bioproject_tallies[bioproject, 'sample_number'] <- length(all_files)
  
  
  raw_in[[bioproject]] <- mclapply(all_files,
                     function(file) {
                       cov_in <- read.table(file, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
                       
                       colnames(cov_in) <- c('genome', 'mean', 'rpkm', 'trimmed_mean', 'relabun', 'covered_bases', 'length')
                       
                       cov_in$breadth <- cov_in$covered_bases / cov_in$length
                       
                       orig_col <- colnames(cov_in)
                       
                       mgs_id <- gsub('\\..*$', '', basename(file))
                       cov_in$mgs_sample <- mgs_id
                       
                       cov_in$bioproject <- bioproject
                       
                       return(cov_in[, c('bioproject', 'mgs_sample', orig_col)])

                      },
                     mc.cores=64)
}

# Realized after running CoverM that the IDs here are actually run IDs, rather than sample IDs.
# So kept only one run per sample (whichever has more genomes present).
derep_runs <- function(run_ids, metadata_tab, coverm_tab, presence_breadth_cutoff) {
  if (length(setdiff(run_ids, rownames(metadata_tab))) > 0) { stop("Error - non-overlapping run IDs!") }
  metadata_tab <- metadata_tab[run_ids, ]
  sample_ids <- metadata_tab[run_ids, "Sample"]
  sample_tallies <- table(sample_ids)
  dup_samples <- names(sample_tallies)[which(sample_tallies > 1)]
  
  nondup_samples <- names(sample_tallies)[which(sample_tallies == 1)]
  
  nodup_runs <- character()
  for (nondup_sample in nondup_samples) {
    tmp <- metadata_tab[which(metadata_tab$Sample == nondup_sample), ]
    if (nrow(tmp) != 1) { stop("Error!")}
    nodup_runs <- c(nodup_runs, rownames(tmp))
  }
  
  rep_dup_runs <- character()
  for (dup_sample in dup_samples) {
    sample_metadata <- metadata_tab[which(metadata_tab$Sample == dup_sample), ]
    
    dup_run_present <- numeric()
    for (dup_run in rownames(sample_metadata)) {
      dup_run_coverm <- coverm_tab[which(coverm_tab$mgs_sample == dup_run), ]
      dup_run_present <- c(dup_run_present, length(which(dup_run_coverm$breadth >= presence_breadth_cutoff)))
    }
    best_sample <- rownames(sample_metadata)[which.max(dup_run_present)]
    rep_dup_runs <- c(rep_dup_runs, best_sample)
  }
  
  run_map <- data.frame(sample_name = c(dup_samples, nondup_samples),
                        run_ids = c(rep_dup_runs, nodup_runs))
  
  if (length(which(duplicated(run_map$sample_name))) > 0) { stop("Error - duplicated output") }
  if (length(which(duplicated(run_map$run_ids))) > 0) { stop("Error - duplicated output") }
  
  rownames(run_map) <- run_map$run_ids
  
  return(run_map)
}

PRJEB1787_coverm <- do.call(rbind, raw_in$PRJEB1787)
PRJEB1787_run_ids <- unique(PRJEB1787_coverm$mgs_sample)
PRJEB1787_metadata <- read.table("/mfs/gdouglas/projects/ocean_mags/metadata/PRJEB1787_metadata.csv",
                                 header = TRUE, sep = ",", row.names = 1, stringsAsFactors = FALSE)
PRJEB1787_sample_run_map <- derep_runs(run_ids = PRJEB1787_run_ids,
                                       metadata_tab = PRJEB1787_metadata,
                                       coverm_tab = PRJEB1787_coverm,
                                       presence_breadth_cutoff = 0.3)
PRJEB1787_coverm <- PRJEB1787_coverm[which(PRJEB1787_coverm$mgs_sample %in% PRJEB1787_sample_run_map$run_ids), ]
PRJEB1787_coverm$mgs_sample_recoded <- PRJEB1787_sample_run_map[PRJEB1787_coverm$mgs_sample, "sample_name"]
PRJEB1787_coverm$mgs_sample <- PRJEB1787_coverm$mgs_sample_recoded

write.table(x = PRJEB1787_sample_run_map,
            file = "/mfs/gdouglas/projects/ocean_mags/coverm/combined_tables/Tara_dedup_id_mapping/PRJEB1787.tsv",
            col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

PRJEB9740_coverm <- do.call(rbind, raw_in$PRJEB9740)
PRJEB9740_run_ids <- unique(PRJEB9740_coverm$mgs_sample)
PRJEB9740_metadata <- read.table("/mfs/gdouglas/projects/ocean_mags/metadata/PRJEB9740_metadata.csv",
                                 header = TRUE, sep = ",", row.names = 1, stringsAsFactors = FALSE)
PRJEB9740_sample_run_map <- derep_runs(run_ids = PRJEB9740_run_ids,
                                       metadata_tab = PRJEB9740_metadata,
                                       coverm_tab = PRJEB9740_coverm,
                                       presence_breadth_cutoff = 0.3)
PRJEB9740_coverm <- PRJEB9740_coverm[which(PRJEB9740_coverm$mgs_sample %in% PRJEB9740_sample_run_map$run_ids), ]
PRJEB9740_coverm$mgs_sample_recoded <- PRJEB9740_sample_run_map[PRJEB9740_coverm$mgs_sample, "sample_name"]
PRJEB9740_coverm$mgs_sample <- PRJEB9740_coverm$mgs_sample_recoded

write.table(x = PRJEB9740_sample_run_map,
            file = "/mfs/gdouglas/projects/ocean_mags/coverm/combined_tables/Tara_dedup_id_mapping/PRJEB9740.tsv",
            col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

PRJEB6608_coverm <- do.call(rbind, raw_in$PRJEB6608)
PRJEB6608_run_ids <- unique(PRJEB6608_coverm$mgs_sample)
PRJEB6608_metadata <- read.table("/mfs/gdouglas/projects/ocean_mags/metadata/PRJEB6608_metadata.csv",
                                 header = TRUE, sep = ",", row.names = 1, stringsAsFactors = FALSE)
PRJEB6608_sample_run_map <- derep_runs(run_ids = PRJEB6608_run_ids,
                                       metadata_tab = PRJEB6608_metadata,
                                       coverm_tab = PRJEB6608_coverm,
                                       presence_breadth_cutoff = 0.15)
PRJEB6608_coverm <- PRJEB6608_coverm[which(PRJEB6608_coverm$mgs_sample %in% PRJEB6608_sample_run_map$run_ids), ]
PRJEB6608_coverm$mgs_sample_recoded <- PRJEB6608_sample_run_map[PRJEB6608_coverm$mgs_sample, "sample_name"]
PRJEB6608_coverm$mgs_sample <- PRJEB6608_coverm$mgs_sample_recoded

write.table(x = PRJEB6608_sample_run_map,
            file = "/mfs/gdouglas/projects/ocean_mags/coverm/combined_tables/Tara_dedup_id_mapping/PRJEB6608.tsv",
            col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

# Split into metaG and metaT datasets.
combined_dna <- rbind(PRJEB1787_coverm,
                      PRJEB9740_coverm)
combined_rna <- PRJEB6608_coverm

# Create presence tables.
dna_presence <- prep_presence_tab(cov_tab=combined_dna, presence_cutoff=0.30)
rna_presence <- prep_presence_tab(cov_tab=combined_rna, presence_cutoff=0.15)

# Create RPKM tables, based on the same samples/genomes as in the presence tables.
dna_rpkm <- reshape2::dcast(data = combined_dna,
                            formula = mgs_sample ~ genome,
                            value.var = "rpkm",
                            fill = 0)

rna_rpkm <- reshape2::dcast(data = combined_rna,
                            formula = mgs_sample ~ genome,
                            value.var = "rpkm",
                            fill = 0)

dna_rpkm <- dna_rpkm[which(dna_rpkm$mgs_sample %in% rownames(dna_presence)), c("mgs_sample", colnames(dna_presence))]
rna_rpkm <- rna_rpkm[which(rna_rpkm$mgs_sample %in% rownames(rna_presence)), c("mgs_sample", colnames(rna_presence))]

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

write_table_w_rowcol(x = dna_presence, outfile = "~/projects/ocean_mags/coverm/combined_tables/metaG_presence.tsv")
write_table_w_rowcol(x = rna_presence, outfile = "~/projects/ocean_mags/coverm/combined_tables/metaT_presence.tsv")

write.table(x = dna_rpkm, file = "~/projects/ocean_mags/coverm/combined_tables/metaG_rpkm.tsv",
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

write.table(x = rna_rpkm, file = "~/projects/ocean_mags/coverm/combined_tables/metaT_rpkm.tsv",
                     sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
