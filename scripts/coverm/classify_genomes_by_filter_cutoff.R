rm(list = ls(all.names = TRUE))

# Get genomes classified by whether they are associated with free-living *or* less-filtered samples based
# both prevalence and relative abundance.

# Also get this split restricted to samples in the "Baltic Sea" and "Hawaii bloom" datasets.
# Note: unfortunately there were very few genomes that were specifically linked to these datasets.

source('~/scripts/ocean_mag_hgt/scripts/coverm/classify_genomes_by_filter_cutoff_function.R')

freeliving_map <- read.table('/mfs/gdouglas/projects/ocean_mags/metadata/OceanDNA_supp_metadata/OceanDNA_metadata_water_freeliving.tsv',
                             header=TRUE, sep='\t', stringsAsFactors = FALSE)

# Re-code one funky sample name to match presence table.
freeliving_map[which(freeliving_map$sample_name == 'ERS492821_ERS492814'), 'sample_name'] <- 'ERS492814'

lessfiltered_map <- read.table('/mfs/gdouglas/projects/ocean_mags/metadata/OceanDNA_supp_metadata/OceanDNA_metadata_water_notparticledepleted.tsv',
                               header=TRUE, sep='\t', stringsAsFactors = FALSE)

classify_genomes_filt_cutoff(freeliving_samples = freeliving_map$sample_name,
                             lessfiltered_samples = lessfiltered_map$sample_name,
                             presence_filepath = '~/projects/ocean_mags/networks/combined_tables/metaG_presence_allsamples.tsv.gz',
                             rpkm_filepath = '~/projects/ocean_mags/networks/combined_tables/metaG_rpkm_allsamples.tsv.gz',
                             freeliv_outfile = '/mfs/gdouglas/projects/ocean_hgt_zenodo/taxa_subsets/Tara_genomes_freeliving_associated.tsv',
                             lessfilt_outfile = '/mfs/gdouglas/projects/ocean_hgt_zenodo/taxa_subsets/Tara_genomes_lessfiltered_associated.tsv',
                             num_cores=64)

# classify_genomes_filt_cutoff(freeliving_samples = freeliving_map[which(freeliving_map$division == 'Baltic Sea'), 'sample_name'],
#                              lessfiltered_samples = lessfiltered_map[which(lessfiltered_map$division == 'Baltic Sea'), 'sample_name'],
#                              presence_filepath = '~/projects/ocean_mags/networks/combined_tables/metaG_presence_allsamples.tsv.gz',
#                              rpkm_filepath = '~/projects/ocean_mags/networks/combined_tables/metaG_rpkm_allsamples.tsv.gz',
#                              freeliv_outfile = '/mfs/gdouglas/projects/ocean_hgt_zenodo/taxa_subsets/Tara_Baltic_Sea_genomes_freeliving_associated.tsv',
#                              lessfilt_outfile = '/mfs/gdouglas/projects/ocean_hgt_zenodo/taxa_subsets/Tara_Baltic_Sea_genomes_lessfiltered_associated.tsv',
#                              num_cores=64)
#
# classify_genomes_filt_cutoff(freeliving_samples = freeliving_map[which(freeliving_map$division == 'Hawaii bloom'), 'sample_name'],
#                              lessfiltered_samples = lessfiltered_map[which(lessfiltered_map$division == 'Hawaii bloom'), 'sample_name'],
#                              presence_filepath = '~/projects/ocean_mags/networks/combined_tables/metaG_presence_allsamples.tsv.gz',
#                              rpkm_filepath = '~/projects/ocean_mags/networks/combined_tables/metaG_rpkm_allsamples.tsv.gz',
#                              freeliv_outfile = '/mfs/gdouglas/projects/ocean_hgt_zenodo/taxa_subsets/Tara_Hawaii_bloom_genomes_freeliving_associated.tsv',
#                              lessfilt_outfile = '/mfs/gdouglas/projects/ocean_hgt_zenodo/taxa_subsets/Tara_Hawaii_bloom_genomes_lessfiltered_associated.tsv',
#                              num_cores=64)
