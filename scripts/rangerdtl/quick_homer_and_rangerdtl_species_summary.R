rm(list = ls(all.names = TRUE))

# Quick summary of HoMer and RANGER-DTL combined summary across all tested species.
library(ggplot2)

analyzed_species <- read.table('/mfs/gdouglas/projects/ocean_mags/water_mag_analysis/species_DTL_analyses/species_to_analyze.txt',
                               stringsAsFactors = FALSE, sep = '', header = FALSE)$V1

num_genomes <- as.integer()
summaries <- list()

for (species in analyzed_species) {
  
  genome_map_filepath <- paste('/mfs/gdouglas/projects/ocean_mags/water_mag_analysis/species_DTL_analyses/homer_prep_map_only/',
                               species,
                               'map_out/genome_ids.tsv',
                               sep = '/')
  
  num_genomes <- c(num_genomes, nrow(read.table(genome_map_filepath, sep = "\t", header = FALSE)))
  
  summary_filepath <- paste('/mfs/gdouglas/projects/ocean_mags/water_mag_analysis/species_DTL_analyses/homer_rangerdtl_summaries',
                            species,
                            'transfers.tsv',
                            sep = '/')
  
  summaries[[species]] <- read.table(summary_filepath, header = TRUE, sep = '\t', stringsAsFactors = FALSE)

}

num_total_gene_transfers <- sapply(summaries, function(x) { nrow(x) })
num_single_gene_transfers <- sapply(summaries, function(x) { nrow(x[which(x$hgt_instance == "single"), ]) })
num_multi_gene_transfer_events <- sapply(summaries, function(x) { length(unique(x[which(x$hgt_instance != "single"), "hgt_instance"])) })
num_multi_gene_transfer_genes <- sapply(summaries, function(x) { length(x[which(x$hgt_instance != "single"), "hgt_instance"]) })

full_summary <- data.frame(species = analyzed_species,
                           num_genomes = num_genomes,
                           num_total_gene_transfers = num_total_gene_transfers,
                           num_single_gene_transfers = num_single_gene_transfers,
                           num_multi_gene_transfer_events = num_multi_gene_transfer_events,
                           num_multi_gene_transfer_genes = num_multi_gene_transfer_genes)

full_summary$percent_multi_genes <- (full_summary$num_multi_gene_transfer_genes / full_summary$num_total_gene_transfers) * 100

full_summary$mean_multi_size <- full_summary$num_multi_gene_transfer_genes / full_summary$num_multi_gene_transfer_events

ggplot(data = full_summary, aes(x = num_genomes, y = num_total_gene_transfers)) +
  geom_point() +
  theme_bw() +
  xlab("Number of genomes") +
  ylab("Number of total gene transfers")

ggplot(data = full_summary, aes(x = num_genomes, y = num_single_gene_transfers)) +
  geom_point() +
  theme_bw() +
  xlab("Number of genomes") +
  ylab("Number of non-HGMT gene transfers")

ggplot(data = full_summary, aes(x = num_genomes, y = percent_multi_genes)) +
  geom_point() +
  theme_bw() +
  xlab("Number of genomes") +
  ylab("Percent transfers in HMGTs")

ggplot(data = full_summary, aes(x = num_genomes, y = num_multi_gene_transfer_events)) +
  geom_point() +
  theme_bw() +
  xlab("Number of genomes") +
  ylab("Number HMGT transfer events")

write.table(x = full_summary,
            file = "/mfs/gdouglas/projects/ocean_mags/water_mag_analysis/species_DTL_analyses/homer_rangerdtl_combined_summary.tsv",
            quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
