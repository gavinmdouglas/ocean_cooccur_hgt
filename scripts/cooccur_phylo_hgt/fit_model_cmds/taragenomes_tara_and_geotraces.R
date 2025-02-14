rm(list = ls(all.names = TRUE))

source('~/scripts/ocean_mag_hgt/scripts/cooccur_phylo_hgt/function_compute_hgt_cooccur_LMs.R')

# Tara genomes, split by Tara and geotraces samples only:

null_out <- compute_hgt_cooccur_LMs(prepped_table_path="/mfs/gdouglas/projects/ocean_mags/networks/dataset_subset/metaG_hyperg_cooccur.geotraces.combined.prepped.tsv.gz",
                                    outfolder="/mfs/gdouglas/projects/ocean_mags/glmm_working/tara_genomes/clusterbased_geotraces_samples_hyperg/",
                                    num_cores=8)

null_out <- compute_hgt_cooccur_LMs(prepped_table_path="/mfs/gdouglas/projects/ocean_mags/networks/dataset_subset/metaG_hyperg_cooccur.Tara.combined.prepped.tsv.gz",
                                    outfolder="/mfs/gdouglas/projects/ocean_mags/glmm_working/tara_genomes/clusterbased_tara_samples_hyperg/",
                                    num_cores=8)
