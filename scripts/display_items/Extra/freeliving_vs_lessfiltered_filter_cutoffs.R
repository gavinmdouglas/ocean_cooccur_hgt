rm(list = ls(all.names = TRUE))

library(ggplot2)
library(cowplot)

presence_tab <- read.table('/mfs/gdouglas/projects/ocean_mags/networks/combined_tables/metaG_presence_allsamples.tsv.gz',
                           header=TRUE, sep='\t', row.names=1)

alldata_meta <- read.table('/mfs/gdouglas/projects/ocean_mags/metadata/OceanDNA_supp_metadata/Supp_File_S1_water_samples_clean.tsv',
                           header=TRUE, sep='\t', stringsAsFactors = FALSE)

freeliving_meta <- read.table('/mfs/gdouglas/projects/ocean_mags/metadata/OceanDNA_supp_metadata/OceanDNA_metadata_water_freeliving.tsv',
                              header=TRUE, sep='\t', stringsAsFactors = FALSE)

lessfiltered_meta <- read.table('/mfs/gdouglas/projects/ocean_mags/metadata/OceanDNA_supp_metadata/OceanDNA_metadata_water_notparticledepleted.tsv',
                                header=TRUE, sep='\t', stringsAsFactors = FALSE)

rownames(alldata_meta) <- alldata_meta$sample_name
rownames(freeliving_meta) <- freeliving_meta$sample_name
rownames(lessfiltered_meta) <- lessfiltered_meta$sample_name

other_ids <- setdiff(rownames(alldata_meta), c(rownames(freeliving_meta), rownames(lessfiltered_meta), 'ERS492821_ERS492814', 'ERS492814'))

other_meta <- alldata_meta[other_ids, ]

other_meta <- other_meta[intersect(rownames(other_meta), rownames(presence_tab)), ]
freeliving_meta <- freeliving_meta[intersect(rownames(freeliving_meta), c(rownames(presence_tab), 'ERS492821_ERS492814')), ]
lessfiltered_meta <- lessfiltered_meta[intersect(rownames(lessfiltered_meta), rownames(presence_tab)), ]

swap_values <- function(vec, orig, new) {
  if (length(orig) != 1) { stop('Orig not length 1') }
  if (length(new) != 1) { stop('New not length 1') }

  if (is.na(orig)) {
    if (length(which(is.na(vec))) > 0) {
      vec[which(is.na(vec))] <- new
    }
  } else {
    if (length(which(vec == orig)) > 0) {
      vec[which(vec == orig)] <- new
    }
  }

  return(vec)
}

compute_filter_cutoff_tallies <- function(intab, category) {

  lower <- as.character(intab$lower_filter)
  higher <- as.character(intab$upper_filter)

  lower <- swap_values(lower, 'NA', '-')
  higher <- swap_values(higher, 'NA', '-')
  higher <- swap_values(higher, '3', '3.0')
  higher <- swap_values(higher, '5', '5.0')

  cutoffs <- paste(lower, higher, sep='|')
  cutoff_tallies <- table(cutoffs)

  lower_num <- gsub('\\|.*$', '', names(cutoff_tallies))
  lower_num <- swap_values(lower_num, '-', '0')
  lower_num <- as.numeric(lower_num)
  upper_num <- gsub('^.*\\|', '', names(cutoff_tallies))
  upper_num <- swap_values(upper_num, '-', '10000')
  upper_num <- as.numeric(upper_num)

  out_df <- data.frame(category = category,
                       cutoffs=names(cutoff_tallies),
                       tallies = as.integer(cutoff_tallies),
                       lower_num = lower_num,
                       upper_num = upper_num)
  out_df <- out_df[order(out_df$lower_num, out_df$upper_num, decreasing = TRUE), ]

  return(out_df)
}

other_breakdown <- compute_filter_cutoff_tallies(other_meta, category='Other')
freeliving_breakdown <- compute_filter_cutoff_tallies(freeliving_meta, category='Free living')
lessfiltered_breakdown <- compute_filter_cutoff_tallies(lessfiltered_meta, category='Less-filtered')

combined_breakdown <- rbind(other_breakdown, freeliving_breakdown)
combined_breakdown <- rbind(combined_breakdown, lessfiltered_breakdown)

combined_breakdown$cutoffs <- factor(combined_breakdown$cutoffs, levels=combined_breakdown$cutoffs)


ggplot(data = combined_breakdown, aes(y=cutoffs, x=tallies)) +
  geom_col() +
  theme_bw() +
  facet_wrap(~category, scales='free') +
  ylab('Lower filter | Upper filter')