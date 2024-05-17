rm(list = ls(all.names = TRUE))

library(reshape2)
library(parallel)
library(data.table)
library(backbone)

hyperg_cooccur_parallel <- function(in_df,
                                    output_tab) {
  
  hyperg_out <- backbone::hyperg(t(in_df))
  
  hyperg_out_pos <- hyperg_out$positive
  hyperg_out_pos[upper.tri(hyperg_out_pos, diag = TRUE)] <- NA
  
  hyperg_out_pos <- data.frame(hyperg_out_pos) 
  hyperg_out_pos$taxon_i <- rownames(hyperg_out_pos)
  hyperg_out_pos_long <- reshape2::melt(data = hyperg_out_pos, id.vars = "taxon_i", value.name = "P", variable.name = "taxon_j")
  hyperg_out_pos_long <- hyperg_out_pos_long[which(! is.na(hyperg_out_pos_long$P)), ]
  hyperg_out_pos_long$BH <- p.adjust(hyperg_out_pos_long$P, 'BH')
  hyperg_out_pos_long <- hyperg_out_pos_long[which(hyperg_out_pos_long$BH < 0.05), ]
  
  write.table(x = hyperg_out_pos_long,
              file = output_tab,
              sep = "\t",
              quote = FALSE,
              col.names = TRUE,
              row.names = FALSE)
}

metaG_presence <- read.table(file = "/mfs/gdouglas/projects/water_mags/coverm/combined_tables/metaG_presence.tsv.gz",
                             header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)

metaT_presence <- read.table(file = "/mfs/gdouglas/projects/water_mags/coverm/combined_tables/metaT_presence.tsv.gz",
                             header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)

metaG_simple_NULL <- hyperg_cooccur_parallel(in_df = metaG_presence,
                                             output_tab = "/mfs/gdouglas/projects/water_mags/coverm/network_working/metaG_hyperg_cooccur.tsv")

metaT_simple_NULL <- hyperg_cooccur_parallel(in_df = metaT_presence,
                                             output_tab = "/mfs/gdouglas/projects/water_mags/coverm/network_working/metaT_hyperg_cooccur.tsv")
