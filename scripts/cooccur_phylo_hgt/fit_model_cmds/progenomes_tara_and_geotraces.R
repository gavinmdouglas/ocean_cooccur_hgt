rm(list = ls(all.names = TRUE))

source('~/scripts/ocean_mag_hgt/scripts/cooccur_phylo_hgt/function_compute_hgt_cooccur_LMs.R')

null_out <- compute_hgt_cooccur_LMs(prepped_table_path="/mfs/gdouglas/projects/ocean_mags/progenomes_analyses/networks/dataset_subset/metaG_hyperg_cooccur.geotraces.combined.prepped.tsv.gz",
                                    outfolder="/mfs/gdouglas/projects/ocean_mags/glmm_working/progenomes/clusterbased_hyperg_geotraces_samples/",
                                    num_cores=8)

null_out <- compute_hgt_cooccur_LMs(prepped_table_path="/mfs/gdouglas/projects/ocean_mags/progenomes_analyses/networks/dataset_subset/metaG_hyperg_cooccur.Tara.combined.prepped.tsv.gz",
                                    outfolder="/mfs/gdouglas/projects/ocean_mags/glmm_working/progenomes/clusterbased_hyperg_tara_samples/",
                                    num_cores=8)

