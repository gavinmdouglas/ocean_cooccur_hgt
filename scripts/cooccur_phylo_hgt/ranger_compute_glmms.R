rm(list = ls(all.names = TRUE))

source('~/scripts/ocean_mag_hgt/scripts/cooccur_phylo_hgt/function_compute_hgt_cooccur_glmms.R')

# null_out <- compute_hgt_cooccur_glmms(combined_file="/mfs/gdouglas/projects/ocean_mags/networks/allsamples/rangerdtl/metaG_hyperg_cooccur.ranger.combined.tsv.gz",
#                                       cooccur_approach='hyperg',
#                                       hgt_tally_col='ranger_hgt_tallies',
#                                       median_diff_file="/mfs/gdouglas/projects/ocean_mags/networks/combined_tables/allsamples_present_metadata_median_pairwise_diff.tsv.gz",
#                                       outfolder="/mfs/gdouglas/projects/ocean_mags/glmm_working/rangerdtl/hyperg",
#                                       num_cores=8,
#                                       keep_lower_levels=TRUE)
#
# null_out <- compute_hgt_cooccur_glmms(combined_file="/mfs/gdouglas/projects/ocean_mags/networks/allsamples/rangerdtl/metaG_simple_cooccur.ranger.combined.tsv.gz",
#                                       cooccur_approach='simple',
#                                       hgt_tally_col='ranger_hgt_tallies',
#                                       median_diff_file="/mfs/gdouglas/projects/ocean_mags/networks/combined_tables/allsamples_present_metadata_median_pairwise_diff.tsv.gz",
#                                       outfolder="/mfs/gdouglas/projects/ocean_mags/glmm_working/rangerdtl/simple",
#                                       num_cores=8,
#                                       keep_lower_levels=TRUE)
#
# null_out <- compute_hgt_cooccur_glmms(combined_file="/mfs/gdouglas/projects/ocean_mags/networks/allsamples/rangerdtl/metaG_propr_rpkm.ranger.combined.tsv.gz",
#                                       cooccur_approach='propr',
#                                       hgt_tally_col='ranger_hgt_tallies',
#                                       median_diff_file="/mfs/gdouglas/projects/ocean_mags/networks/combined_tables/allsamples_present_metadata_median_pairwise_diff.tsv.gz",
#                                       outfolder="/mfs/gdouglas/projects/ocean_mags/glmm_working/rangerdtl/propr",
#                                       num_cores=8,
#                                       keep_lower_levels=TRUE)

null_out <- compute_hgt_cooccur_glmms(combined_file="/mfs/gdouglas/projects/ocean_mags/networks/filtersplit/rangerdtl/metaG_hyperg_cooccur.freeliv.combined.tsv.gz",
                                      cooccur_approach='hyperg',
                                      hgt_tally_col='ranger_hgt_tallies',
                                      median_diff_file="/mfs/gdouglas/projects/ocean_mags/networks/combined_tables/freeliv_present_metadata_median_pairwise_diff.tsv.gz",
                                      outfolder="/mfs/gdouglas/projects/ocean_mags/glmm_working/rangerdtl/hyperg",
                                      num_cores=8,
                                      suffix='.freeliv.rds',
                                      keep_lower_levels=TRUE)

null_out <- compute_hgt_cooccur_glmms(combined_file="/mfs/gdouglas/projects/ocean_mags/networks/filtersplit/rangerdtl/metaG_simple_cooccur.freeliv.combined.tsv.gz",
                                      cooccur_approach='simple',
                                      hgt_tally_col='ranger_hgt_tallies',
                                      median_diff_file="/mfs/gdouglas/projects/ocean_mags/networks/combined_tables/freeliv_present_metadata_median_pairwise_diff.tsv.gz",
                                      outfolder="/mfs/gdouglas/projects/ocean_mags/glmm_working/rangerdtl/simple",
                                      num_cores=8,
                                      suffix='.freeliv.rds',
                                      keep_lower_levels=TRUE)


null_out <- compute_hgt_cooccur_glmms(combined_file="/mfs/gdouglas/projects/ocean_mags/networks/filtersplit/rangerdtl/metaG_propr_rpkm.freeliv.combined.tsv.gz",
                                      cooccur_approach='propr',
                                      hgt_tally_col='ranger_hgt_tallies',
                                      median_diff_file="/mfs/gdouglas/projects/ocean_mags/networks/combined_tables/freeliv_present_metadata_median_pairwise_diff.tsv.gz",
                                      outfolder="/mfs/gdouglas/projects/ocean_mags/glmm_working/rangerdtl/propr",
                                      num_cores=8,
                                      suffix='.freeliv.rds',
                                      keep_lower_levels=TRUE)

# Filter-split - lessfiltered
null_out <- compute_hgt_cooccur_glmms(combined_file="/mfs/gdouglas/projects/ocean_mags/networks/filtersplit/rangerdtl/metaG_hyperg_cooccur.lessfiltered.combined.tsv.gz",
                                      cooccur_approach='hyperg',
                                      hgt_tally_col='ranger_hgt_tallies',
                                      median_diff_file="/mfs/gdouglas/projects/ocean_mags/networks/combined_tables/lessfiltered_present_metadata_median_pairwise_diff.tsv.gz",
                                      outfolder="/mfs/gdouglas/projects/ocean_mags/glmm_working/rangerdtl/hyperg",
                                      num_cores=8,
                                      suffix='.lessfiltered.rds',
                                      keep_lower_levels=TRUE)

null_out <- compute_hgt_cooccur_glmms(combined_file="/mfs/gdouglas/projects/ocean_mags/networks/filtersplit/rangerdtl/metaG_simple_cooccur.lessfiltered.combined.tsv.gz",
                                      cooccur_approach='simple',
                                      hgt_tally_col='ranger_hgt_tallies',
                                      median_diff_file="/mfs/gdouglas/projects/ocean_mags/networks/combined_tables/lessfiltered_present_metadata_median_pairwise_diff.tsv.gz",
                                      outfolder="/mfs/gdouglas/projects/ocean_mags/glmm_working/rangerdtl/simple",
                                      num_cores=8,
                                      suffix='.lessfiltered.rds',
                                      keep_lower_levels=TRUE)

null_out <- compute_hgt_cooccur_glmms(combined_file="/mfs/gdouglas/projects/ocean_mags/networks/filtersplit/rangerdtl/metaG_propr_rpkm.lessfiltered.combined.tsv.gz",
                                      cooccur_approach='propr',
                                      hgt_tally_col='ranger_hgt_tallies',
                                      median_diff_file="/mfs/gdouglas/projects/ocean_mags/networks/combined_tables/lessfiltered_present_metadata_median_pairwise_diff.tsv.gz",
                                      outfolder="/mfs/gdouglas/projects/ocean_mags/glmm_working/rangerdtl/propr",
                                      num_cores=8,
                                      suffix='.lessfiltered.rds',
                                      keep_lower_levels=TRUE)
