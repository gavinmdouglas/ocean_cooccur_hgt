library(NetCoMi)
library(rlang)
library(reshape2)
library(parallel)
library(data.table)
library(backbone)
library(Matrix)
library(reshape2)
library(parallel)
library(data.table)

hyperg_cooccur_parallel <- function(in_df,
                                    output_gzip_file,
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

  gzfile_connection <-  gzfile(output_gzip_file, "w")

  write.table(x = hyperg_out_pos_long,
              file = gzfile_connection,
              sep = "\t",
              quote = FALSE,
              col.names = TRUE,
              row.names = FALSE)

  close(gzfile_connection)

}


run_netcomi_propr <- function(in_tab, random_seed, outfile_gizpped) {
  output <- netConstruct(in_tab,
                         measure = "propr",
                         sparsMethod="threshold",
                         thresh = 0,
                         zeroMethod="multRepl",
                         verbose = 3,
                         seed = random_seed)

  colnames(output$edgelist1)[1] <- "taxon_i"
  colnames(output$edgelist1)[2] <- "taxon_j"

  gzfile_connection <- gzfile(outfile_gizpped, "w")

  write.table(x = output$edgelist1,
              file = gzfile_connection,
              sep = "\t",
              quote = FALSE,
              col.names = TRUE,
              row.names = FALSE)

  close(gzfile_connection)
}


simple_cooccur_parallel <- function(in_df,
                                    ncores,
                                    tmp_dir,
                                    output_gzip_file) {

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

  gzfile_connection <-  gzfile(output_gzip_file, "w")

  write.table(x = simple_cooccur_combined,
              file = gzfile_connection,
              sep = "\t",
              quote = FALSE,
              col.names = TRUE,
              row.names = FALSE)

  close(gzfile_connection)
}

