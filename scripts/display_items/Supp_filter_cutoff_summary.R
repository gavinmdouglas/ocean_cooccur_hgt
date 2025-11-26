rm(list = ls(all.names = TRUE))

# Summarize number of samples by their filter cut-off groupings for the three general groupings:
# (1) small, (2) less-filtered, and (3) large

library(cowplot)
library(ggplot2)

info <- read.table("/mfs/gdouglas/projects/ocean_mags/metadata/OceanDNA_supp_metadata/Supp_File_S1_water_samples_clean.tsv",
                   header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Then for filters.
filters <- data.frame(lower=info$lower_filter, upper=info$upper_filter)

filters$Category <- "Other"
filters$Category[which(filters$upper <= 3)] <- "Free-living"
filters$Category[which(filters$upper >= 20 | is.na(filters$upper))] <- "Free-living and particle-attached"

filters$lower <- as.character(filters$lower)
filters$upper <- as.character(filters$upper)

filters$upper[which(filters$upper == '3')] <- '3.0'
filters$upper[which(filters$upper == '5')] <- '5.0'

filters$lower[is.na(filters$lower)] <- 'NA'
filters$upper[is.na(filters$upper)] <- 'NA'

filters$pair <- paste(filters$lower, filters$upper, sep = ' - ')

pair_tallies <- table(filters$pair)
filters <- filters[-which(duplicated(filters)), ]

filters$upper[which(filters$upper == '-')] <- '100000'
filters$lower[which(filters$lower == '-')] <- '0'

filters$upper <- as.numeric(filters$upper)
filters$lower <- as.numeric(filters$lower)

filters$Tally <- as.numeric(pair_tallies[filters$pair])

filters <- filters[order(filters$upper, filters$lower, decreasing = FALSE), ]

filters$pair <- factor(filters$pair, levels = filters$pair)

filters <- filters[-which(filters$Category == "Other"), ]

filters$Category <- factor(filters$Category, levels=c("Free-living", "Free-living and particle-attached"))

oceanDNA_filtercutoffs <- ggplot(data=filters, aes(x=Tally, y=pair)) +
                              geom_col() +
                              theme_bw() +
                              geom_text(aes(label = Tally), hjust = -0.2, size = 3) +
                              facet_wrap(~ Category, scales = "free") +
                              xlab('Metagenomics sample count') +
                              ylab('Lower filter | Upper filter') +
                              scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
                              theme(text = element_text(colour="black"),
                                    axis.text = element_text(colour="black"),
                                    strip.background = element_blank(),
                                    strip.text = element_text(colour = "black"))


# And then for upper fraction samples.

env_data <- read.table("/mfs/gdouglas/projects/ocean_hgt_zenodo/mapfiles/Tara_PANGAEA_env_data.tsv.gz",
                       header=TRUE, sep = '\t', stringsAsFactors = FALSE, row.names = 1)

env_data <- env_data[, -which(colnames(env_data) %in% c("Sample_ID_BioSamples_accession_number", "Sample_ID_ENA_sample_accession_number"))]

env_data_upper <- env_data[which(env_data$Fraction_lower_µm_used_on_board_to_prepare_samp >= 3 & env_data$Fraction_upper_µm_used_on_board_to_prepare_samp == 20), ]

# Only plot for the 43 samples that were analyzed:
combined_summary <- read.table('/mfs/gdouglas/projects/ocean_hgt_zenodo/hgt_prev_analyses/Tara_fraction_sample_matched_HGT_prev.tsv.gz',
                               header=TRUE, sep = '\t', stringsAsFactors = FALSE)

combined_summary <- combined_summary[which(combined_summary$lower_total_genomes >= 10 & combined_summary$upper_total_genomes >= 10), ]

env_data_upper <- env_data_upper[rownames(env_data_upper) %in% combined_summary$upper_sample, ]

upper <- data.frame(lower=env_data_upper$Fraction_lower_µm_used_on_board_to_prepare_samp,
                    upper=env_data_upper$Fraction_upper_µm_used_on_board_to_prepare_samp)

upper$lower <- as.character(upper$lower)
upper$upper <- as.character(upper$upper)

upper$lower[which(upper$lower == '3')] <- '3.0'
upper$lower[which(upper$lower == '5')] <- '5.0'

upper$pair <- paste(upper$lower, upper$upper, sep = '|')

upper_pair_tallies <- table(upper$pair)
upper <- upper[-which(duplicated(upper)), ]

upper$Tally <- as.numeric(upper_pair_tallies[upper$pair])

upper <- upper[order(upper$upper, upper$lower, decreasing = FALSE), ]

upper$pair <- factor(upper$pair, levels = upper$pair)

upper_filtercutoffs <- ggplot(data=upper, aes(x=Tally, y=pair)) +
  geom_col() +
  theme_bw() +
  geom_text(aes(label = Tally), hjust = -0.2, size = 3) +
  xlab('Metagenomics sample count') +
  ylab('Lower filter | Upper filter') +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  theme(text = element_text(colour="black"),
        axis.text = element_text(colour="black"))

lower_row <- plot_grid(NULL, upper_filtercutoffs, NULL, rel_widths = c(1, 2, 1), ncol = 3)

combined_plot <- plot_grid(oceanDNA_filtercutoffs, lower_row, nrow=2, rel_heights = c(3.5, 1))

combined_plot <- plot_grid(combined_plot, labels=c('b'))

ggsave(plot = combined_plot,
       filename = "/mfs/gdouglas/scripts/ocean_cooccur_hgt/display_items/Supp_filter_cutoff_breakdown_RAW.pdf",
       device = "pdf", width = 8, height = 5, units = "in", dpi=600)
