
identify_enriched_categories <- function(genes,
                                         background,
                                         gene_to_category_map,
                                         min_num = 0,
                                         to_ignore = c()) {

  enrichments_out <- data.frame(matrix(NA, nrow = length(gene_to_category_map), ncol = 8))
  rownames(enrichments_out) <- names(gene_to_category_map)
  colnames(enrichments_out) <- c("category", "genes_num_category", "genes_num_other",
                                 "background_num_category", "background_num_other", "OR", "p", "fdr")

  enrichments_out[names(gene_to_category_map), "category"] <- names(gene_to_category_map)

  for (category in rownames(enrichments_out)) {

    if (category %in% to_ignore) { next }

    genes_num_category <- length(which(genes %in% gene_to_category_map[[category]]))
    genes_num_other <- length(genes) - genes_num_category

    background_num_category <- length(which(background %in% gene_to_category_map[[category]]))
    background_num_other <- length(background) - background_num_category

    count_table <- matrix(c(genes_num_category, genes_num_other, background_num_category, background_num_other), nrow = 2, ncol = 2)

    if (min(count_table) < min_num) { next }

    fisher_out <- fisher.test(count_table)

    enrichments_out[category, c("genes_num_category",
                                "genes_num_other",
                                "background_num_category",
                                "background_num_other", "p")] <- c(genes_num_category,
                                                                   genes_num_other,
                                                                   background_num_category,
                                                                   background_num_other,
                                                                   fisher_out$p.value)
    if (genes_num_other > 0) {
      ratio_numer <- genes_num_category / genes_num_other
    } else {
      ratio_numer <- genes_num_category / 1
    }

    if (background_num_other == 0) {
      ratio_denom <- 1
    } else if(background_num_category == 0) {
      ratio_denom <- 1 / background_num_other
    } else {
      ratio_denom <- background_num_category / background_num_other
    }

    enrichments_out[category, "OR"] <- ratio_numer / ratio_denom
  }

  if (length(which(rowSums(is.na(enrichments_out)) > 1)) > 0) {
    enrichments_out <- enrichments_out[-which(rowSums(is.na(enrichments_out)) > 1), ]
  }

  enrichments_out$fdr <- p.adjust(enrichments_out$p, "fdr")

  rownames(enrichments_out) <- NULL

  return(enrichments_out)

}
