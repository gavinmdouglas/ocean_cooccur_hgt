rm(list = ls(all.names = TRUE))

library(reshape2)
library(parallel)
library(data.table)

simple_cooccur_parallel <- function(in_df,
                                    ncores,
                                    tmp_dir,
                                    output_tab) {
  
  presence_vecs <- lapply(in_df, function(x) { which(x == 1) })
  
  presence_vecs_lengths <- sapply(presence_vecs, length)
  
  # Get index combinations to loop over.
  i_combinations <- combn(1:length(presence_vecs), 2)

  i_combinations_list <- mclapply(1:ncol(i_combinations),
                                  function(i) { as.integer(i_combinations[, i]) },
                                  mc.cores = ncores)
  
  total_comb <- length(i_combinations_list)
  
  combo_batch_ranges <- c(seq(from=1, to=length(i_combinations_list), by=1000000), length(i_combinations_list))
  
  obs_comb <- 0
  for (batch_i in 1:(length(combo_batch_ranges) - 1)) {

    start_index <- combo_batch_ranges[batch_i]
    end_index <- combo_batch_ranges[batch_i + 1] - 1
    
    if (batch_i == (length(combo_batch_ranges) - 1)) {
      end_index <- combo_batch_ranges[batch_i + 1]
    }
    
    simple_cooccur_raw <- mclapply(i_combinations_list[start_index:end_index],
                                   function(combo) {
                                     i = combo[1]
                                     j = combo[2]
                                     num_intersect = length(intersect(presence_vecs[[i]],
                                                                      presence_vecs[[j]]))
                                     
                                     min_num = min(presence_vecs_lengths[[i]],
                                                   presence_vecs_lengths[[j]])
                                     
                                     return(list(taxoni=i,
                                                 taxonj=j,
                                                 num_intersect=num_intersect,
                                                 min_num=min_num,
                                                 simple_cooccur=num_intersect/min_num))
                                   },
                                   mc.cores = ncores)
    
    simple_cooccur <- data.frame(data.table::rbindlist(simple_cooccur_raw))
    
    obs_comb <- obs_comb + nrow(simple_cooccur)
    
    simple_cooccur <- simple_cooccur[which(simple_cooccur$simple_cooccur > 0), ]
    simple_cooccur$taxoni <- names(presence_vecs)[simple_cooccur$taxoni]
    simple_cooccur$taxonj <- names(presence_vecs)[simple_cooccur$taxonj]
    
    saveRDS(object = simple_cooccur,
            file = paste(tmp_dir, "/simple_cooccur_set", as.character(batch_i), ".rds", sep = ""))
  }
  
  message("Expected ", total_comb, " comparisons")
  message("Observed ", obs_comb, " comparisons")
  if (obs_comb != total_comb) {
    stop("Mismatch in combo comparisons.")
  }

  output_rds_files <- list.files(tmp_dir, full.names = TRUE)
  
  simple_cooccur_combined <- readRDS(output_rds_files[1])
  
  for (i in 2:length(output_rds_files)) {
    simple_cooccur_combined <- rbind(simple_cooccur_combined, readRDS(output_rds_files[i]))
  }
  
  write.table(x = simple_cooccur_combined,
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

metaG_simple_NULL <- simple_cooccur_parallel(in_df = metaG_presence,
                                             ncores = 64,
                                             tmp_dir = "/mfs/gdouglas/projects/water_mags/coverm/network_working/metaG_simple_cooccur_tmp/",
                                             output_tab = "/mfs/gdouglas/projects/water_mags/coverm/network_working/metaG_simple_cooccur.tsv")

metaT_simple_NULL <- simple_cooccur_parallel(in_df = metaT_presence,
                                             ncores = 64,
                                             tmp_dir = "/mfs/gdouglas/projects/water_mags/coverm/network_working/metaT_simple_cooccur_tmp/",
                                             output_tab = "/mfs/gdouglas/projects/water_mags/coverm/network_working/metaT_simple_cooccur.tsv")
    