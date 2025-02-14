compute_mean_median_sample_vals <- function(info_tab, presence_tab, outprefix) {

  vars_of_interest <- c("lower_filter", "upper_filter", "depth", "latitude",
                        "longitude", "temperature", "oxygen", "salinity")

  mean_tab <- data.frame(matrix(NA,
                                nrow = ncol(presence_tab),
                                ncol = length(vars_of_interest)))
  colnames(mean_tab) <- vars_of_interest
  rownames(mean_tab) <- colnames(presence_tab)

  for (taxon in colnames(presence_tab)) {
    samples_w_taxon <- rownames(presence_tab)[which(presence_tab[, taxon] > 0)]
    if (length(samples_w_taxon) < 10) { stop("Error - should be at least 10 MGS samples!") }

    info_taxon_subset <- info_tab[samples_w_taxon, ]

    for (var_of_interest in vars_of_interest) {
      vec <- info_taxon_subset[, var_of_interest]
      if (length(which(! is.na(vec))) > 0) {
        mean_tab[taxon, var_of_interest] <- mean(vec, na.rm = TRUE)
      }
    }
  }

  mean_outfile <- paste(outprefix, 'mean_by_sample.tsv.gz', sep = '_')
  mean_gzfile_connection <-  gzfile(mean_outfile, "w")
  write.table(x = mean_tab,
              file = mean_gzfile_connection,
              col.names = NA, row.names = TRUE, sep = "\t", quote = FALSE)
  close(mean_gzfile_connection)

  median_tab <- data.frame(matrix(NA,
                                  nrow = ncol(presence_tab),
                                  ncol = length(vars_of_interest)))
  colnames(median_tab) <- vars_of_interest
  rownames(median_tab) <- colnames(presence_tab)

  for (taxon in colnames(presence_tab)) {
    samples_w_taxon <- rownames(presence_tab)[which(presence_tab[, taxon] > 0)]
    if (length(samples_w_taxon) < 10) { stop("Error - should be at least 10 MGS samples!") }

    info_taxon_subset <- info_tab[samples_w_taxon, ]

    for (var_of_interest in vars_of_interest) {
      vec <- info_taxon_subset[, var_of_interest]
      if (length(which(! is.na(vec))) > 0) {
        median_tab[taxon, var_of_interest] <- median(vec, na.rm = TRUE)
      }
    }
  }

  median_outfile <- paste(outprefix, 'median_by_sample.tsv.gz', sep = '_')
  median_gzfile_connection <- gzfile(median_outfile, "w")
  write.table(x = median_tab,
              file = median_gzfile_connection,
              col.names = NA, row.names = TRUE, sep = "\t", quote = FALSE)
  close(median_gzfile_connection)
}
