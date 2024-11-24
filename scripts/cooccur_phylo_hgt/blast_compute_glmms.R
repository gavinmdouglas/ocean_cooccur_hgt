rm(list = ls(all.names = TRUE))

source('~/scripts/ocean_mag_hgt/scripts/cooccur_phylo_hgt/function_compute_hgt_cooccur_glmms.R')

# null_out <- compute_hgt_cooccur_glmms(combined_file="/mfs/gdouglas/projects/ocean_mags/networks/allsamples/blast/metaG_hyperg_cooccur.allsamples.combined.tsv.gz",
#                                       cooccur_approach='hyperg',
#                                       hgt_tally_col='both_hit_count',
#                                       median_diff_file="/mfs/gdouglas/projects/ocean_mags/networks/combined_tables/allsamples_present_metadata_median_pairwise_diff.tsv.gz",
#                                       outfolder="/mfs/gdouglas/projects/ocean_mags/glmm_working/blast/hyperg",
#                                       num_cores=8)
#
# null_out <- compute_hgt_cooccur_glmms(combined_file="/mfs/gdouglas/projects/ocean_mags/networks/allsamples/blast/metaG_simple_cooccur.allsamples.combined.tsv.gz",
#                                       cooccur_approach='simple',
#                                       hgt_tally_col='both_hit_count',
#                                       median_diff_file="/mfs/gdouglas/projects/ocean_mags/networks/combined_tables/allsamples_present_metadata_median_pairwise_diff.tsv.gz",
#                                       outfolder="/mfs/gdouglas/projects/ocean_mags/glmm_working/blast/simple",
#                                       num_cores=8)
#
# null_out <- compute_hgt_cooccur_glmms(combined_file="/mfs/gdouglas/projects/ocean_mags/networks/allsamples/blast/metaG_propr_rpkm.allsamples.combined.tsv.gz",
#                                       cooccur_approach='propr',
#                                       hgt_tally_col='both_hit_count',
#                                       median_diff_file="/mfs/gdouglas/projects/ocean_mags/networks/combined_tables/allsamples_present_metadata_median_pairwise_diff.tsv.gz",
#                                       outfolder="/mfs/gdouglas/projects/ocean_mags/glmm_working/blast/propr",
#                                       num_cores=8)


# Filter-split - freeliv
null_out <- compute_hgt_cooccur_glmms(combined_file="/mfs/gdouglas/projects/ocean_mags/networks/filtersplit/blast/metaG_hyperg_cooccur.freeliv.combined.tsv.gz",
                                      cooccur_approach='hyperg',
                                      hgt_tally_col='both_hit_count',
                                      median_diff_file="/mfs/gdouglas/projects/ocean_mags/networks/combined_tables/freeliv_present_metadata_median_pairwise_diff.tsv.gz",
                                      outfolder="/mfs/gdouglas/projects/ocean_mags/glmm_working/blast/hyperg",
                                      num_cores=8,
                                      suffix='.freeliv.rds')

null_out <- compute_hgt_cooccur_glmms(combined_file="/mfs/gdouglas/projects/ocean_mags/networks/filtersplit/blast/metaG_simple_cooccur.freeliv.combined.tsv.gz",
                                      cooccur_approach='simple',
                                      hgt_tally_col='both_hit_count',
                                      median_diff_file="/mfs/gdouglas/projects/ocean_mags/networks/combined_tables/freeliv_present_metadata_median_pairwise_diff.tsv.gz",
                                      outfolder="/mfs/gdouglas/projects/ocean_mags/glmm_working/blast/simple",
                                      num_cores=8,
                                      suffix='.freeliv.rds')


null_out <- compute_hgt_cooccur_glmms(combined_file="/mfs/gdouglas/projects/ocean_mags/networks/filtersplit/blast/metaG_propr_rpkm.freeliv.combined.tsv.gz",
                                      cooccur_approach='propr',
                                      hgt_tally_col='both_hit_count',
                                      median_diff_file="/mfs/gdouglas/projects/ocean_mags/networks/combined_tables/freeliv_present_metadata_median_pairwise_diff.tsv.gz",
                                      outfolder="/mfs/gdouglas/projects/ocean_mags/glmm_working/blast/propr",
                                      num_cores=8,
                                      suffix='.freeliv.rds')

# Filter-split - lessfiltered
null_out <- compute_hgt_cooccur_glmms(combined_file="/mfs/gdouglas/projects/ocean_mags/networks/filtersplit/blast/metaG_hyperg_cooccur.lessfiltered.combined.tsv.gz",
                                      cooccur_approach='hyperg',
                                      hgt_tally_col='both_hit_count',
                                      median_diff_file="/mfs/gdouglas/projects/ocean_mags/networks/combined_tables/lessfiltered_present_metadata_median_pairwise_diff.tsv.gz",
                                      outfolder="/mfs/gdouglas/projects/ocean_mags/glmm_working/blast/hyperg",
                                      num_cores=8,
                                      suffix='.lessfiltered.rds')

null_out <- compute_hgt_cooccur_glmms(combined_file="/mfs/gdouglas/projects/ocean_mags/networks/filtersplit/blast/metaG_simple_cooccur.lessfiltered.combined.tsv.gz",
                                      cooccur_approach='simple',
                                      hgt_tally_col='both_hit_count',
                                      median_diff_file="/mfs/gdouglas/projects/ocean_mags/networks/combined_tables/lessfiltered_present_metadata_median_pairwise_diff.tsv.gz",
                                      outfolder="/mfs/gdouglas/projects/ocean_mags/glmm_working/blast/simple",
                                      num_cores=8,
                                      suffix='.lessfiltered.rds')

null_out <- compute_hgt_cooccur_glmms(combined_file="/mfs/gdouglas/projects/ocean_mags/networks/filtersplit/blast/metaG_propr_rpkm.lessfiltered.combined.tsv.gz",
                                      cooccur_approach='propr',
                                      hgt_tally_col='both_hit_count',
                                      median_diff_file="/mfs/gdouglas/projects/ocean_mags/networks/combined_tables/lessfiltered_present_metadata_median_pairwise_diff.tsv.gz",
                                      outfolder="/mfs/gdouglas/projects/ocean_mags/glmm_working/blast/propr",
                                      num_cores=8,
                                      suffix='.lessfiltered.rds')

