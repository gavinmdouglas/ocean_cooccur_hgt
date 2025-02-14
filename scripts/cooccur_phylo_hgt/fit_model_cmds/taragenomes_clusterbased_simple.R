rm(list = ls(all.names = TRUE))

source('~/scripts/ocean_mag_hgt/scripts/cooccur_phylo_hgt/function_compute_hgt_cooccur_LMs.R')

null_out <- compute_hgt_cooccur_LMs(prepped_table_path="/mfs/gdouglas/projects/ocean_mags/networks/allsamples/clusterbased/metaG_simple_cooccur.allsamples.combined.prepped.tsv.gz",
                                    outfolder="/mfs/gdouglas/projects/ocean_mags/glmm_working/tara_genomes/clusterbased_allsamples/simple/",
                                    num_cores=8)

null_out <- compute_hgt_cooccur_LMs_w_filtergroups(prepped_table_path="/mfs/gdouglas/projects/ocean_mags/networks/allsamples/clusterbased/metaG_simple_cooccur.allsamples.combined.prepped.tsv.gz",
                                                   outfolder="/mfs/gdouglas/projects/ocean_mags/glmm_working/tara_genomes/clusterbased_allsamples/simple/",
                                                   num_cores=8,
                                                   freeliving_genomes_file='/mfs/gdouglas/projects/ocean_hgt_zenodo/taxa_subsets/Tara_genomes_freeliving_associated.tsv.gz',
                                                   lessfiltered_genomes_file='/mfs/gdouglas/projects/ocean_hgt_zenodo/taxa_subsets/Tara_genomes_lessfiltered_associated.tsv.gz')
