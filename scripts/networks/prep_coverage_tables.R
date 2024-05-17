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

bioprojects <- c('PRJEB1787', 'PRJEB97740', 'PRJEB6608')

bioproject_tallies <- data.frame(bioprojects = bioprojects,
                                 sample_number = NA)
rownames(bioproject_tallies) <- bioprojects

raw_in <- list()

for (bioproject in bioprojects) {
  
  all_files <- list.files(paste0("~/projects/water_mags/coverm/", bioproject), 
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

# Split into metaG and metaT datasets.
combined_dna <- rbind(do.call(rbind, raw_in$PRJEB1787),
                      do.call(rbind, raw_in$PRJEB97740))
combined_rna <- do.call(rbind, raw_in$PRJEB6608)

# Create presence tables.
dna_presence <- prep_presence_tab(cov_tab=combined_dna, presence_cutoff=0.30)

rna_presence <- prep_presence_tab(cov_tab=combined_rna, presence_cutoff=0.15)
rna_presence <- rna_presence[which(rowSums(rna_presence) >= 30), ]

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

write_table_w_rowcol(x = dna_presence, outfile = "~/projects/water_mags/coverm/combined_tables/metaG_presence.tsv")
write_table_w_rowcol(x = rna_presence, outfile = "~/projects/water_mags/coverm/combined_tables/metaT_presence.tsv")

write.table(x = dna_rpkm, file = "~/projects/water_mags/coverm/combined_tables/metaG_rpkm.tsv",
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

write.table(x = rna_rpkm, file = "~/projects/water_mags/coverm/combined_tables/metaT_rpkm.tsv",
                     sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
