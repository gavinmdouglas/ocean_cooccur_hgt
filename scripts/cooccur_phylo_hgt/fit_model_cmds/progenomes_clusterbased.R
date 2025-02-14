rm(list = ls(all.names = TRUE))

source('~/scripts/ocean_mag_hgt/scripts/cooccur_phylo_hgt/function_compute_hgt_cooccur_LMs.R')

# Progenomes analysis -- all samples.

null_out <- compute_hgt_cooccur_LMs(prepped_table_path="/mfs/gdouglas/projects/ocean_mags/progenomes_analyses/networks/allsamples/clusterbased/metaG_hyperg_cooccur.allsamples.combined.prepped.tsv.gz",
                                    outfolder="/mfs/gdouglas/projects/ocean_mags/glmm_working/progenomes/clusterbased_hyperg_allsamples/",
                                    num_cores=8)

null_out <- compute_hgt_cooccur_LMs_w_filtergroups(prepped_table_path="/mfs/gdouglas/projects/ocean_mags/progenomes_analyses/networks/allsamples/clusterbased/metaG_hyperg_cooccur.allsamples.combined.prepped.tsv.gz",
                                                   outfolder="/mfs/gdouglas/projects/ocean_mags/glmm_working/progenomes/clusterbased_hyperg_allsamples/",
                                                   num_cores=8,
                                                   freeliving_genomes_file='/mfs/gdouglas/projects/ocean_hgt_zenodo/taxa_subsets/progenomes_freeliving_associated.tsv.gz',
                                                   lessfiltered_genomes_file='/mfs/gdouglas/projects/ocean_hgt_zenodo/taxa_subsets/progenomes_lessfiltered_associated.tsv.gz')
