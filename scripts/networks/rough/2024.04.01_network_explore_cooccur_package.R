rm(list = ls(all.names = TRUE))

library(reshape2)
library(parallel)
library(data.table)
library(cooccur)

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

# Then analyze data - starting with the DNA projects.
combined_dna$bioproject <- factor(combined_dna$bioproject, levels = bioprojects)

combined_dna$present <- 0

combined_dna[which(combined_dna$breadth >= 0.25), "present"] <- 1

dna_presence_matrix <- reshape2::dcast(data = combined_dna,
                                      formula = mgs_sample ~ genome,
                                      value.var = "present",
                                      fill = 0)

rownames(dna_presence_matrix) <- dna_presence_matrix$mgs_sample
dna_presence_matrix <- dna_presence_matrix[, which(colnames(dna_presence_matrix) != "mgs_sample")]

dna_sample_tallies <- rowSums(dna_presence_matrix)
dna_presence_matrix_filt <- dna_presence_matrix[which(dna_sample_tallies >= 50), ]

dna_taxon_presence_tallies <- colSums(dna_presence_matrix_filt)
dna_presence_matrix_filt <- dna_presence_matrix_filt[, which(dna_taxon_presence_tallies >= 5)]

dna_presence_matrix_filt_t <- t(dna_presence_matrix_filt)

backbone_out <- backbone::hyperg(dna_presence_matrix_filt_t)

backbone_out_pos <- backbone_out$positive
backbone_out_pos[upper.tri(backbone_out_pos, diag = TRUE)] <- NA

backbone_out_pos_prep <- data.frame(backbone_out_pos) 
backbone_out_pos_prep$taxon_i <- rownames(backbone_out_pos_prep)
backbone_out_pos_long <- reshape2::melt(data = backbone_out_pos_prep, id.vars = "taxon_i", value.name = "P", variable.name = "taxon_j")
backbone_out_pos_long <- backbone_out_pos_long[which(! is.na(backbone_out_pos_long$P)), ]
backbone_out_pos_long$BH <- p.adjust(backbone_out_pos_long$P, 'BH')
backbone_out_pos_long_sig <- backbone_out_pos_long[which(backbone_out_pos_long$BH < 0.05), ]

saveRDS(object = backbone_out_pos_long_sig,
        file = "/mfs/gdouglas/projects/water_mags/coverm/network_working/dna_backbone_cooccur.rds")
