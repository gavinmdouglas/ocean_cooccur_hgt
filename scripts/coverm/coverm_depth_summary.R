rm(list = ls(all.names = TRUE))

library(parallel)

all_files <- c(list.files("/mfs/gdouglas/projects/ocean_mags/coverm/additional_output", pattern = "cov.gz", full.names = TRUE),
               list.files("/mfs/gdouglas/projects/ocean_mags/coverm/additional_OceanDNA_round2", pattern = "cov.gz", full.names = TRUE))

sample_ids <- sapply(all_files, function(x) { gsub('.TARA_MAGs.cov.gz', '', basename(x)) })
names(sample_ids) <- NULL
length(which(duplicated(sample_ids)))

num_reads_mapped <- unlist(mclapply(all_files,
                                    function(x) {
                                      in_tab <- read.table(x, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
                                      in_tab <- in_tab[which(in_tab$Genome != 'unmapped'), ]
                                      read_count_col_i <- grep('Read.Count', colnames(in_tab))
                                      if (length(read_count_col_i) != 1) {
                                        stop('Should only find one col!')
                                      }
                                      sum(in_tab[, read_count_col_i], na.rm = TRUE)
                                    }, mc.cores=64))

coverm_reads_mapped <- data.frame(mgs_sample=sample_ids,
                                  num_reads_mapped=num_reads_mapped)

called_present <- read.table('~/projects/ocean_mags/networks/combined_tables/metaG_presence_allsamples.tsv.gz',
                             header=TRUE, sep = '\t', stringsAsFactors = FALSE, row.names = 1)

excluded_samples <- setdiff(coverm_reads_mapped$mgs_sample, rownames(called_present))

coverm_reads_mapped[which(coverm_reads_mapped$mgs_sample %in% excluded_samples), ]

coverm_reads_mapped <- coverm_reads_mapped[which(! coverm_reads_mapped$mgs_sample %in% excluded_samples), ]

summary(coverm_reads_mapped$num_reads_mapped)
sd(coverm_reads_mapped$num_reads_mapped)

