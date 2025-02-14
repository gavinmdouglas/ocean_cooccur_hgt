rm(list = ls(all.names = TRUE))

source('~/scripts/ocean_mag_hgt/scripts/cooccur_phylo_hgt/function_compute_hgt_cooccur_glmms.R')

null_out <- compute_hgt_cooccur_w_subgroups(combined_file="/mfs/gdouglas/projects/ocean_mags/networks/allsamples/blast/metaG_hyperg_cooccur.allsamples.combined.tsv.gz",
                                            cooccur_approach='hyperg',
                                            hgt_tally_col='both_hit_count',
                                            median_diff_file="/mfs/gdouglas/projects/ocean_mags/networks/combined_tables/allsamples_present_metadata_median_pairwise_diff.tsv.gz",
                                            outfolder="/mfs/gdouglas/projects/ocean_mags/glmm_working/blast/hyperg",
                                            num_cores=8)

