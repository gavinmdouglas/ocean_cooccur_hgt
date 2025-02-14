rm(list = ls(all.names = TRUE))

# Get tallies of how many pairwise inter-taxa level comparisons there were for each level.
# Do for MAGs subset for MetaG and MetaT tables separately.

num_pairwise_comparisons_cross_highest_category <- function(tax_table, decreasing_tax_columns) {
  
  level_tallies <- rep(0, length(decreasing_tax_columns))
  names(level_tallies) <- decreasing_tax_columns
  
  for (i in 1:(nrow(tax_table) - 1)) {
    
    tax_table_remainder <- tax_table[(i + 1):nrow(tax_table), , drop = FALSE]
    
    for (tax_level in decreasing_tax_columns) {
      
      diff_index <- which(tax_table[i, tax_level] != tax_table_remainder[, tax_level])
      
      # Skip if no differences.
      if (length(diff_index) == 0) {
        next
      } else {
        level_tallies[tax_level] <- level_tallies[tax_level] + length(diff_index)
        tax_table_remainder <- tax_table_remainder[-diff_index, , drop = FALSE]
      }
      
    }
    
  }
  
  return(level_tallies)
  
}

taxa_levels <- c("Strain", "Species", "Genus", "Family", "Order", "Class", "Phylum", "Domain")

taxa <- read.table("/mfs/gdouglas/projects/water_mags/water_mag_analysis/mapfiles/MAG_taxa_breakdown.tsv.gz",
                   header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = 2)

metaG_presence <- read.table("/mfs/gdouglas/projects/water_mags/coverm/combined_tables/metaG_presence.tsv.gz",
                             header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = 1)

metaT_presence <- read.table("/mfs/gdouglas/projects/water_mags/coverm/combined_tables/metaT_presence.tsv.gz",
                             header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = 1)

taxa_metaG <- taxa[colnames(metaG_presence), ]
taxa_metaT <- taxa[colnames(metaT_presence), ]

taxa_metaG_num_comparisons <- num_pairwise_comparisons_cross_highest_category(tax_table = taxa_metaG,
                                                                              decreasing_tax_columns = rev(taxa_levels))

taxa_metaT_num_comparisons <- num_pairwise_comparisons_cross_highest_category(tax_table = taxa_metaT,
                                                                              decreasing_tax_columns = rev(taxa_levels))

taxa_metaG_num_comparisons_df <- data.frame(Level = names(taxa_metaG_num_comparisons), Num_comparisons = taxa_metaG_num_comparisons)

taxa_metaT_num_comparisons_df <- data.frame(Level = names(taxa_metaT_num_comparisons), Num_comparisons = taxa_metaT_num_comparisons)

write.table(x = taxa_metaG_num_comparisons_df,
            file = "/mfs/gdouglas/projects/water_mags/num_comparisons/metaG_taxa_compare_per_inter.level.tsv",
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

write.table(x = taxa_metaT_num_comparisons_df,
            file = "/mfs/gdouglas/projects/water_mags/num_comparisons/metaT_taxa_compare_per_inter.level.tsv",
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
