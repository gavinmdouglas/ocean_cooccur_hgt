classify_genomes_filt_cutoff <- function(freeliving_samples,
                                         lessfiltered_samples,
                                         presence_filepath,
                                         rpkm_filepath,
                                         freeliv_outfile,
                                         lessfilt_outfile,
                                         num_cores=64) {

  freeliving_samples = freeliving_map$sample_name
  lessfiltered_samples = lessfiltered_map$sample_name
  presence_filepath = '~/projects/ocean_mags/networks/combined_tables/metaG_presence_allsamples.tsv.gz'
  rpkm_filepath = '~/projects/ocean_mags/networks/combined_tables/metaG_rpkm_allsamples.tsv.gz'
  freeliv_outfile = '/mfs/gdouglas/projects/ocean_hgt_zenodo/taxa_subsets/Tara_genomes_freeliving_associated.tsv'
  lessfilt_outfile = '/mfs/gdouglas/projects/ocean_hgt_zenodo/taxa_subsets/Tara_genomes_lessfiltered_associated.tsv'
  num_cores=64

  presence <- read.table(presence_filepath, header=TRUE, sep='\t', stringsAsFactors = FALSE, row.names=1)

  rpkm <- read.table(rpkm_filepath, header=TRUE, sep='\t', stringsAsFactors = FALSE, row.names=1)

  freeliving_samples <- intersect(freeliving_samples, rownames(presence))
  lessfiltered_samples <- intersect(lessfiltered_samples, rownames(presence))

  num_freeliving_genomes_per_sample <- rowSums(presence[freeliving_samples, ])
  num_lessfiltered_genomes_per_sample <- rowSums(presence[lessfiltered_samples, ])

  prop_freeliving <- length(freeliving_samples) / nrow(presence)
  prop_lessfiltered <- length(lessfiltered_samples) / nrow(presence)

  count_and_rpkm_breakdown <- function(genome) {
    genome_samples <- rownames(presence)[which(presence[, genome] > 0)]
    num_samples <- length(genome_samples)

    freeliving_genome_samples <- intersect(genome_samples, freeliving_samples)
    lessfiltered_genome_samples <- intersect(genome_samples, lessfiltered_samples)
    other_genome_samples <- setdiff(genome_samples, c(freeliving_genome_samples, lessfiltered_genome_samples))

    freeliving_median_rpkm <- median(rpkm[freeliving_genome_samples, genome])
    lessfiltered_median_rpkm <- median(rpkm[lessfiltered_genome_samples, genome])
    other_median_rpkm <- median(rpkm[other_genome_samples, genome])

    highest_median_rpkm_group <- c('freeliv', 'lessfilt', 'other')[which.max(c(freeliving_median_rpkm, lessfiltered_median_rpkm, other_median_rpkm))]

    return(data.frame(genome = genome,
                      num_samples = num_samples,
                      freeliving_num_samples = length(freeliving_genome_samples),
                      lessfiltered_num_samples = length(lessfiltered_genome_samples),
                      freeliving_median_rpkm = freeliving_median_rpkm,
                      lessfiltered_median_rpkm = lessfiltered_median_rpkm,
                      other_median_rpkm = other_median_rpkm,
                      highest_median_rpkm_group = highest_median_rpkm_group))
  }

  taxa_count_breakdown_raw <- parallel::mclapply(colnames(presence), count_and_rpkm_breakdown, mc.cores=num_cores)

  taxa_count_breakdown <- do.call(rbind, taxa_count_breakdown_raw)

  taxa_count_breakdown$other_samples <- taxa_count_breakdown$num_samples - (taxa_count_breakdown$freeliving_num_samples + taxa_count_breakdown$lessfiltered_num_samples)

  taxa_count_breakdown$freeliving_genome_prop <- taxa_count_breakdown$freeliving_num_samples / taxa_count_breakdown$num_samples
  taxa_count_breakdown$lessfiltered_genome_prop <- taxa_count_breakdown$lessfiltered_num_samples / taxa_count_breakdown$num_samples

  genomes_mainly_freeliving <- taxa_count_breakdown$genome[which(taxa_count_breakdown$freeliving_genome_prop > 0.75 & taxa_count_breakdown$highest_median_rpkm_group == 'freeliv')]
  genomes_mainly_lessfiltered <- taxa_count_breakdown$genome[which(taxa_count_breakdown$lessfiltered_genome_prop > 0.75 & taxa_count_breakdown$highest_median_rpkm_group == 'lessfilt')]

  write.table(x = genomes_mainly_freeliving, file = freeliv_outfile, col.names = FALSE, row.names = FALSE, quote = FALSE)
  write.table(x = genomes_mainly_lessfiltered, file = lessfilt_outfile, col.names = FALSE, row.names = FALSE, quote = FALSE)

}
