rm(list = ls(all.names = TRUE))

# Prepare columns of input tables, to make them fast to re-run.
prep_input_tab <- function(combined_file,
                           cooccur_approach,
                           hgt_tally_col,
                           median_diff_file,
                           outfile,
                           keep_lower_levels) {

  if (! cooccur_approach %in% c('hyperg', 'simple', 'propr')) { stop('Co-occur approach must be hyperg, simple, or propr.') }
  if (length(grep(cooccur_approach, combined_file)) == 0) { stop('Co-occur approach not present in combined_file?') }
  if (!grepl("\\.gz$", outfile)) { stop("Error: outfile must have a .gz extension") }

  print(combined_file)

  combined_info <- read.table(combined_file, header=TRUE, sep="\t", stringsAsFactors = FALSE, row.names = 1)

  if (! keep_lower_levels) {
    combined_info <- combined_info[which(! combined_info$diff_tax_level %in% c('Species', 'Strain')), ]
  }
  combined_info <- combined_info[which(rowSums(is.na(combined_info)) == 0), ]

  if (length(grep('.ranger.', combined_file)) == 0) {
    combined_info <- combined_info[which(combined_info$diff_tax_level %in% c('Kingdom', 'Domain', 'Phylum', 'Class', 'Order', 'Genus', 'Species', 'Strain')), ]
  }

  if (cooccur_approach == 'hyperg') {
    combined_info$cooccur_ratio <- combined_info$cooccur_obs / combined_info$cooccur_exp
    combined_info$cooccur <- 0
    combined_info$cooccur[which(combined_info$cooccur_BH < 0.05 & combined_info$cooccur_ratio > 1)] <- 1
  } else if (cooccur_approach == 'simple') {
    combined_info <- combined_info[which(! is.na(combined_info$cooccur_simple_cooccur)), ]
    combined_info$cooccur <- bestNormalize::orderNorm(x = combined_info$cooccur_simple_cooccur)$x.t
  } else if (cooccur_approach == 'propr') {
    combined_info <- combined_info[which(! is.na(combined_info$cooccur_asso)), ]
    combined_info$cooccur <- bestNormalize::orderNorm(x = combined_info$cooccur_asso)$x.t
  }

  combined_info$hgt <- 0
  combined_info$hgt[which(combined_info[, hgt_tally_col] > 0)] <- 1

  combined_info$tip_dist_orderedNorm <- bestNormalize::orderNorm(x = combined_info$tip_dist)$x.t

  # Add in metadata as well.
  pairwise_metadata_dist <- read.table(median_diff_file, header=TRUE, sep="\t", stringsAsFactors = FALSE, row.names = 1)

  combined_info$depth_orderedNorm <-  bestNormalize::orderNorm(x = pairwise_metadata_dist[rownames(combined_info), 'depth'])$x.t
  combined_info$latitude_orderedNorm <-  bestNormalize::orderNorm(x = pairwise_metadata_dist[rownames(combined_info), 'latitude'])$x.t
  combined_info$longitude_orderedNorm <-  bestNormalize::orderNorm(x = pairwise_metadata_dist[rownames(combined_info), 'longitude'])$x.t
  combined_info$temperature_orderedNorm <-  bestNormalize::orderNorm(x = pairwise_metadata_dist[rownames(combined_info), 'temperature'])$x.t
  combined_info$oxygen_orderedNorm <-  bestNormalize::orderNorm(x = pairwise_metadata_dist[rownames(combined_info), 'oxygen'])$x.t
  combined_info$salinity_orderedNorm <-  bestNormalize::orderNorm(x = pairwise_metadata_dist[rownames(combined_info), 'salinity'])$x.t

  combined_info <- combined_info[which(rowSums(is.na(combined_info)) == 0), ]

  write.table(x = combined_info, file = gzfile(outfile), sep = '\t', quote=FALSE, col.names = TRUE, row.names = FALSE)

}

# Prepare columns of input tables, to make them fast to re-run.

# # For Tara genomes, all-samples analysis, with cluster-based HGT.
# null_out <- prep_input_tab(combined_file="/mfs/gdouglas/projects/ocean_mags/networks/allsamples/clusterbased/metaG_hyperg_cooccur.allsamples.combined.tsv.gz",
#                                       cooccur_approach='hyperg',
#                                       hgt_tally_col='both_gene_count',
#                                       median_diff_file="/mfs/gdouglas/projects/ocean_mags/networks/combined_tables/allsamples_present_metadata_median_pairwise_diff.tsv.gz",
#                                       outfile="/mfs/gdouglas/projects/ocean_mags/networks/allsamples/clusterbased/metaG_hyperg_cooccur.allsamples.combined.prepped.tsv.gz",
#                                       keep_lower_levels = FALSE)
#
# null_out <- prep_input_tab(combined_file="/mfs/gdouglas/projects/ocean_mags/networks/allsamples/clusterbased/metaG_simple_cooccur.allsamples.combined.tsv.gz",
#                                       cooccur_approach='simple',
#                                       hgt_tally_col='both_gene_count',
#                                       median_diff_file="/mfs/gdouglas/projects/ocean_mags/networks/combined_tables/allsamples_present_metadata_median_pairwise_diff.tsv.gz",
#                                       outfile="/mfs/gdouglas/projects/ocean_mags/networks/allsamples/clusterbased/metaG_simple_cooccur.allsamples.combined.prepped.tsv.gz",
#                                       keep_lower_levels = FALSE)

# null_out <- prep_input_tab(combined_file="/mfs/gdouglas/projects/ocean_mags/networks/allsamples/clusterbased/metaG_propr_rpkm.allsamples.combined.tsv.gz",
#                                       cooccur_approach='propr',
#                                       hgt_tally_col='both_gene_count',
#                                       median_diff_file="/mfs/gdouglas/projects/ocean_mags/networks/combined_tables/allsamples_present_metadata_median_pairwise_diff.tsv.gz",
#                                       outfile="/mfs/gdouglas/projects/ocean_mags/networks/allsamples/clusterbased/metaG_propr_rpkm.allsamples.combined.prepped.tsv.gz",
#                                       keep_lower_levels = FALSE)

# As above, but for BLAST output.
# null_out <- prep_input_tab(combined_file="/mfs/gdouglas/projects/ocean_mags/networks/allsamples/blast/metaG_hyperg_cooccur.allsamples.combined.tsv.gz",
#                            cooccur_approach='hyperg',
#                            hgt_tally_col='both_hit_count',
#                            median_diff_file="/mfs/gdouglas/projects/ocean_mags/networks/combined_tables/allsamples_present_metadata_median_pairwise_diff.tsv.gz",
#                            outfile="/mfs/gdouglas/projects/ocean_mags/networks/allsamples/blast/metaG_hyperg_cooccur.allsamples.combined.prepped.tsv.gz",
#                            keep_lower_levels = FALSE)
#
# null_out <- prep_input_tab(combined_file="/mfs/gdouglas/projects/ocean_mags/networks/allsamples/blast/metaG_simple_cooccur.allsamples.combined.tsv.gz",
#                            cooccur_approach='simple',
#                            hgt_tally_col='both_hit_count',
#                            median_diff_file="/mfs/gdouglas/projects/ocean_mags/networks/combined_tables/allsamples_present_metadata_median_pairwise_diff.tsv.gz",
#                            outfile="/mfs/gdouglas/projects/ocean_mags/networks/allsamples/blast/metaG_simple_cooccur.allsamples.combined.prepped.tsv.gz",
#                           keep_lower_levels = FALSE)

# null_out <- prep_input_tab(combined_file="/mfs/gdouglas/projects/ocean_mags/networks/allsamples/blast/metaG_propr_rpkm.allsamples.combined.tsv.gz",
#                            cooccur_approach='propr',
#                            hgt_tally_col='both_hit_count',
#                            median_diff_file="/mfs/gdouglas/projects/ocean_mags/networks/combined_tables/allsamples_present_metadata_median_pairwise_diff.tsv.gz",
#                            outfile="/mfs/gdouglas/projects/ocean_mags/networks/allsamples/blast/metaG_propr_rpkm.allsamples.combined.prepped.tsv.gz",
#                            keep_lower_levels = FALSE)

# As above, but for RANGER output.
# null_out <- prep_input_tab(combined_file="/mfs/gdouglas/projects/ocean_mags/networks/allsamples/rangerdtl/metaG_hyperg_cooccur.ranger.combined.tsv.gz",
#                            cooccur_approach='hyperg',
#                            hgt_tally_col='ranger_hgt_tallies',
#                            median_diff_file="/mfs/gdouglas/projects/ocean_mags/networks/combined_tables/allsamples_present_metadata_median_pairwise_diff.tsv.gz",
#                            outfile="/mfs/gdouglas/projects/ocean_mags/networks/allsamples/rangerdtl/metaG_hyperg_cooccur.allsamples.combined.prepped.tsv.gz",
#                            keep_lower_levels = TRUE)
#
# null_out <- prep_input_tab(combined_file="/mfs/gdouglas/projects/ocean_mags/networks/allsamples/rangerdtl/metaG_simple_cooccur.ranger.combined.tsv.gz",
#                            cooccur_approach='simple',
#                            hgt_tally_col='ranger_hgt_tallies',
#                            median_diff_file="/mfs/gdouglas/projects/ocean_mags/networks/combined_tables/allsamples_present_metadata_median_pairwise_diff.tsv.gz",
#                            outfile="/mfs/gdouglas/projects/ocean_mags/networks/allsamples/rangerdtl/metaG_simple_cooccur.allsamples.combined.prepped.tsv.gz",
#                            keep_lower_levels = TRUE)
#
# null_out <- prep_input_tab(combined_file="/mfs/gdouglas/projects/ocean_mags/networks/allsamples/rangerdtl/metaG_propr_rpkm.ranger.combined.tsv.gz",
#                            cooccur_approach='propr',
#                            hgt_tally_col='ranger_hgt_tallies',
#                            median_diff_file="/mfs/gdouglas/projects/ocean_mags/networks/combined_tables/allsamples_present_metadata_median_pairwise_diff.tsv.gz",
#                            outfile="/mfs/gdouglas/projects/ocean_mags/networks/allsamples/rangerdtl/metaG_propr_rpkm.allsamples.combined.prepped.tsv.gz",
#                            keep_lower_levels = TRUE)
#
# # Tara genomes, split by Tara and geotraces samples only:
# null_out <- prep_input_tab(combined_file="/mfs/gdouglas/projects/ocean_mags/networks/dataset_subset/metaG_hyperg_cooccur.geotraces.combined.tsv.gz",
#                            cooccur_approach='hyperg',
#                            hgt_tally_col='both_gene_count',
#                            median_diff_file="/mfs/gdouglas/projects/ocean_mags/networks/combined_tables/geotraces_present_metadata_median_pairwise_diff.tsv.gz",
#                            outfile="/mfs/gdouglas/projects/ocean_mags/networks/dataset_subset/metaG_hyperg_cooccur.geotraces.combined.prepped.tsv.gz",
#                            keep_lower_levels = FALSE)
#
# null_out <- prep_input_tab(combined_file="/mfs/gdouglas/projects/ocean_mags/networks/dataset_subset/metaG_hyperg_cooccur.Tara.combined.tsv.gz",
#                            cooccur_approach='hyperg',
#                            hgt_tally_col='both_gene_count',
#                            median_diff_file="/mfs/gdouglas/projects/ocean_mags/networks/combined_tables/Tara_present_metadata_median_pairwise_diff.tsv.gz",
#                            outfile="/mfs/gdouglas/projects/ocean_mags/networks/dataset_subset/metaG_hyperg_cooccur.Tara.combined.prepped.tsv.gz",
#                            keep_lower_levels = FALSE)

# For progenomes analysis, all samples:
null_out <- prep_input_tab(combined_file="/mfs/gdouglas/projects/ocean_mags/progenomes_analyses/networks/allsamples/clusterbased/metaG_hyperg_cooccur.allsamples.combined.tsv.gz",
                           cooccur_approach='hyperg',
                           hgt_tally_col='both_gene_count',
                           median_diff_file="/mfs/gdouglas/projects/ocean_mags/progenomes_analyses/networks/combined_tables/allsamples_present_metadata_median_pairwise_diff.tsv.gz",
                           outfile="/mfs/gdouglas/projects/ocean_mags/progenomes_analyses/networks/allsamples/clusterbased/metaG_hyperg_cooccur.allsamples.combined.prepped.tsv.gz",
                           keep_lower_levels = FALSE)
#
# # For progenomes analysis, split into Tara and geotraces samples only.
# null_out <- prep_input_tab(combined_file="/mfs/gdouglas/projects/ocean_mags/progenomes_analyses/networks/dataset_subset/metaG_hyperg_cooccur.geotraces.combined.tsv.gz",
#                            cooccur_approach='hyperg',
#                            hgt_tally_col='both_gene_count',
#                            median_diff_file="/mfs/gdouglas/projects/ocean_mags/progenomes_analyses/networks/combined_tables/geotraces_present_metadata_median_pairwise_diff.tsv.gz",
#                            outfile="/mfs/gdouglas/projects/ocean_mags/progenomes_analyses/networks/dataset_subset/metaG_hyperg_cooccur.geotraces.combined.prepped.tsv.gz",
#                            keep_lower_levels = FALSE)
#
# null_out <- prep_input_tab(combined_file="/mfs/gdouglas/projects/ocean_mags/progenomes_analyses/networks/dataset_subset/metaG_hyperg_cooccur.Tara.combined.tsv.gz",
#                cooccur_approach='hyperg',
#                hgt_tally_col='both_gene_count',
#                median_diff_file="/mfs/gdouglas/projects/ocean_mags/progenomes_analyses/networks/combined_tables/Tara_present_metadata_median_pairwise_diff.tsv.gz",
#                outfile="/mfs/gdouglas/projects/ocean_mags/progenomes_analyses/networks/dataset_subset/metaG_hyperg_cooccur.Tara.combined.prepped.tsv.gz",
#                keep_lower_levels = FALSE)
#
