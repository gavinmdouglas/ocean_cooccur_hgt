rm(list = ls(all.names = TRUE))

freeliving_map <- read.table('/mfs/gdouglas/projects/ocean_mags/metadata/OceanDNA_supp_metadata/OceanDNA_metadata_water_freeliving.tsv',
                             header=TRUE, sep='\t', stringsAsFactors = FALSE)

lessfiltered_map <- read.table('/mfs/gdouglas/projects/ocean_mags/metadata/OceanDNA_supp_metadata/OceanDNA_metadata_water_notparticledepleted.tsv',
                               header=TRUE, sep='\t', stringsAsFactors = FALSE)

presence <- read.table('~/projects/ocean_mags/networks/combined_tables/metaG_presence_allsamples.tsv.gz',
                       header=TRUE, sep='\t', stringsAsFactors = FALSE, row.names=1)

freeliving_samples <- intersect(freeliving_map$sample_name, rownames(presence))
lessfiltered_samples <- intersect(lessfiltered_map$sample_name, rownames(presence))

num_freeliving_genomes_per_sample <- rowSums(presence[freeliving_samples, ])
num_lessfiltered_genomes_per_sample <- rowSums(presence[lessfiltered_samples, ])

prop_freeliving <- length(freeliving_samples) / nrow(presence)
prop_lessfiltered <- length(lessfiltered_samples) / nrow(presence)

# Make sure that the number of genomes per sample is in same ball-park (on average) for both sample groupings:
boxplot(num_freeliving_genomes_per_sample, num_lessfiltered_genomes_per_sample)

count_breakdown <- function(genome) {
  genome_samples <- rownames(presence)[which(presence[, genome] > 0)]
  num_samples <- length(genome_samples)
  freeliving_num_samples <- length(intersect(genome_samples, freeliving_map$sample_name))
  lessfiltered_num_samples <- length(intersect(genome_samples, lessfiltered_map$sample_name))
  data.frame(genome = genome, num_samples = num_samples, freeliving_num_samples = freeliving_num_samples, lessfiltered_num_samples = lessfiltered_num_samples)
}

taxa_count_breakdown_raw <- parallel::mclapply(colnames(presence), count_breakdown, mc.cores=64)

taxa_count_breakdown <- do.call(rbind, taxa_count_breakdown_raw)

taxa_count_breakdown$other_samples <- taxa_count_breakdown$num_samples - (taxa_count_breakdown$freeliving_num_samples + taxa_count_breakdown$lessfiltered_num_samples)

taxa_count_breakdown$freeliving_genome_prop <- taxa_count_breakdown$freeliving_num_samples / taxa_count_breakdown$num_samples
taxa_count_breakdown$lessfiltered_genome_prop <- taxa_count_breakdown$lessfiltered_num_samples / taxa_count_breakdown$num_samples

genomes_mainly_freeliving <- taxa_count_breakdown$genome[which(taxa_count_breakdown$freeliving_genome_prop > 0.75)]
genomes_mainly_lessfiltered <- taxa_count_breakdown$genome[which(taxa_count_breakdown$lessfiltered_genome_prop > 0.75)]

write.table(x = genomes_mainly_freeliving, file = '/mfs/gdouglas/projects/ocean_hgt_zenodo/cooccur/genomes_freeliving_associated.tsv',
            col.names = FALSE, row.names = FALSE, quote = FALSE)

write.table(x = genomes_mainly_lessfiltered, file = '/mfs/gdouglas/projects/ocean_hgt_zenodo/cooccur/genomes_lessfiltered_associated.tsv',
            col.names = FALSE, row.names = FALSE, quote = FALSE)

