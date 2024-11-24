rm(list = ls(all.names = TRUE))

source('~/scripts/ocean_mag_hgt/scripts/cooccur_phylo_hgt/function_compute_hgt_cooccur_glmms.R')

# null_out <- compute_multimember_hgt_cooccur_glmms(combined_file="/mfs/gdouglas/projects/ocean_mags/networks/allsamples/rangerdtl/metaG_hyperg_cooccur.ranger.combined.tsv.gz",
#                                                   cooccur_approach='hyperg',
#                                                   hgt_tally_col='ranger_hgt_tallies',
#                                                   median_diff_file="/mfs/gdouglas/projects/ocean_mags/networks/combined_tables/allsamples_present_metadata_median_pairwise_diff.tsv.gz",
#                                                   outprefix="/mfs/gdouglas/projects/ocean_mags/glmm_working/multimember/out_rangerdtl_hyperg",
#                                                   keep_lower_levels=TRUE,
#                                                   suffix='.rds')

null_out <- compute_multimember_hgt_cooccur_glmms(combined_file="/mfs/gdouglas/projects/ocean_mags/networks/allsamples/clusterbased/metaG_hyperg_cooccur.allsamples.combined.tsv.gz",
                                                  cooccur_approach='hyperg',
                                                  hgt_tally_col='both_gene_count',
                                                  median_diff_file="/mfs/gdouglas/projects/ocean_mags/networks/combined_tables/allsamples_present_metadata_median_pairwise_diff.tsv.gz",
                                                  outprefix="/mfs/gdouglas/projects/ocean_mags/glmm_working/multimember/out_clusterbased_hyperg")

null_out <- compute_multimember_hgt_cooccur_glmms(combined_file="/mfs/gdouglas/projects/ocean_mags/networks/allsamples/blast/metaG_hyperg_cooccur.allsamples.combined.tsv.gz",
                                                  cooccur_approach='hyperg',
                                                  hgt_tally_col='both_hit_count',
                                                  median_diff_file="/mfs/gdouglas/projects/ocean_mags/networks/combined_tables/allsamples_present_metadata_median_pairwise_diff.tsv.gz",
                                                  outprefix="/mfs/gdouglas/projects/ocean_mags/glmm_working/multimember/out_blast_hyperg")

null_out <- compute_multimember_hgt_cooccur_glmms(combined_file="/mfs/gdouglas/projects/ocean_mags/networks/allsamples/clusterbased/metaG_simple_cooccur.allsamples.combined.tsv.gz",
                                                  cooccur_approach='simple',
                                                  hgt_tally_col='both_gene_count',
                                                  median_diff_file="/mfs/gdouglas/projects/ocean_mags/networks/combined_tables/allsamples_present_metadata_median_pairwise_diff.tsv.gz",
                                                  outprefix="/mfs/gdouglas/projects/ocean_mags/glmm_working/multimember/out_clusterbased_simple")

null_out <- compute_multimember_hgt_cooccur_glmms(combined_file="/mfs/gdouglas/projects/ocean_mags/networks/allsamples/clusterbased/metaG_propr_rpkm.allsamples.combined.tsv.gz",
                                                  cooccur_approach='propr',
                                                  hgt_tally_col='both_gene_count',
                                                  median_diff_file="/mfs/gdouglas/projects/ocean_mags/networks/combined_tables/allsamples_present_metadata_median_pairwise_diff.tsv.gz",
                                                  outprefix="/mfs/gdouglas/projects/ocean_mags/glmm_working/multimember/out_clusterbased_propr")

null_out <- compute_multimember_hgt_cooccur_glmms(combined_file="/mfs/gdouglas/projects/ocean_mags/networks/filtersplit/clusterbased/metaG_hyperg_cooccur.freeliv.combined.tsv.gz",
                                                  cooccur_approach='hyperg',
                                                  hgt_tally_col='both_gene_count',
                                                  median_diff_file="/mfs/gdouglas/projects/ocean_mags/networks/combined_tables/freeliv_present_metadata_median_pairwise_diff.tsv.gz",
                                                  outprefix="/mfs/gdouglas/projects/ocean_mags/glmm_working/multimember/out_clusterbased_hyperg_freeliv")

null_out <- compute_multimember_hgt_cooccur_glmms(combined_file="/mfs/gdouglas/projects/ocean_mags/networks/filtersplit/clusterbased/metaG_hyperg_cooccur.lessfiltered.combined.tsv.gz",
                                                  cooccur_approach='hyperg',
                                                  hgt_tally_col='both_gene_count',
                                                  median_diff_file="/mfs/gdouglas/projects/ocean_mags/networks/combined_tables/lessfiltered_present_metadata_median_pairwise_diff.tsv.gz",
                                                  outprefix="/mfs/gdouglas/projects/ocean_mags/glmm_working/multimember/out_clusterbased_hyperg_lessfiltered")
