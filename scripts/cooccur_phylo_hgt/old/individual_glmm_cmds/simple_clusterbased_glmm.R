rm(list = ls(all.names = TRUE))

source('~/scripts/ocean_mag_hgt/scripts/cooccur_phylo_hgt/function_compute_hgt_cooccur_glmms.R')

# Increased itermax and evalmax, because the model didn't converge.
null_out <- compute_hgt_cooccur_w_subgroups(combined_file="/mfs/gdouglas/projects/ocean_mags/networks/allsamples/clusterbased/metaG_simple_cooccur.allsamples.combined.tsv.gz",
                                            cooccur_approach='simple',
                                            hgt_tally_col='both_gene_count',
                                            median_diff_file="/mfs/gdouglas/projects/ocean_mags/networks/combined_tables/allsamples_present_metadata_median_pairwise_diff.tsv.gz",
                                            outfolder="/mfs/gdouglas/projects/ocean_mags/glmm_working/clusterbased/simple",
                                            num_cores=8,
                                            itermax=1000,
                                            evalmax=1000)
