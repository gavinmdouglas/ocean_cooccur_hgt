rm(list = ls(all.names = TRUE))

library(reshape2)
library(parallel)
library(data.table)
library(backbone)
library(Matrix)

hyperg_cooccur_parallel <- function(in_df,
                                    output_tab,
                                    min_P = 1,
                                    num_cores = 1) {
  
  hyperg_out <- backbone::hyperg(t(in_df))
  
  hyperg_out_pos <- hyperg_out$positive
  hyperg_out_pos[upper.tri(hyperg_out_pos, diag = TRUE)] <- NA
  
  hyperg_out_pos <- data.frame(hyperg_out_pos) 
  hyperg_out_pos$taxon_i <- rownames(hyperg_out_pos)
  hyperg_out_pos_long <- reshape2::melt(data = hyperg_out_pos, id.vars = "taxon_i", value.name = "P", variable.name = "taxon_j")
  hyperg_out_pos_long <- hyperg_out_pos_long[which(! is.na(hyperg_out_pos_long$P)), ]
  hyperg_out_pos_long$BH <- p.adjust(hyperg_out_pos_long$P, 'BH')
  hyperg_out_pos_long <- hyperg_out_pos_long[which(hyperg_out_pos_long$BH <= min_P), ]
  
  num_comparisons <- nrow(in_df)
  
  process_row <- function(i) {
    taxon1 <- as.character(hyperg_out_pos_long[i, "taxon_i"])
    taxon2 <- as.character(hyperg_out_pos_long[i, "taxon_j"])
    
    set1 <- as.numeric(in_df[, taxon1])
    set2 <- as.numeric(in_df[, taxon2])
    set1_present <- which(set1 > 0)
    set2_present <- which(set2 > 0)
    
    obs <- length(intersect(set1_present, set2_present))
    exp <- ((length(set1_present) / num_comparisons) * (length(set2_present) / num_comparisons)) * num_comparisons
    
    list(exp = exp, obs = obs, ratio = (obs + 1) / (exp + 1))
  }
  
  processed_rows_raw <- mclapply(X = 1:nrow(hyperg_out_pos_long),
                                 FUN = process_row,
                                 mc.cores = num_cores)
  
  hyperg_out_pos_long$exp <- sapply(processed_rows_raw, function(x) { return(x$exp) })
  hyperg_out_pos_long$obs <- sapply(processed_rows_raw, function(x) { return(x$obs) })
  hyperg_out_pos_long$ratio <- sapply(processed_rows_raw, function(x) { return(x$ratio) })
  
  write.table(x = hyperg_out_pos_long,
              file = output_tab,
              sep = "\t",
              quote = FALSE,
              col.names = TRUE,
              row.names = FALSE)
  
}

metaG_presence <- read.table(file = "/mfs/gdouglas/projects/ocean_mags/coverm/combined_tables/metaG_presence_allsamples.tsv.gz",
                             header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)

metaG_hyperg_NULL <- hyperg_cooccur_parallel(in_df = metaG_presence,
                                             output_tab = "/mfs/gdouglas/projects/ocean_mags/coverm/network_working_allsamples/metaG_hyperg_cooccur.allsamples.tsv",
                                             num_cores = 64)
