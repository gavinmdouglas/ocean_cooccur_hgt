rm(list = ls(all.names = TRUE))

# Quick sanity check that the less-filtered bacteria do indeed more often undergo HGT (and that this isn't some artifact in GLMM).

library(ggplot2)
library(reshape2)
library(ComplexHeatmap)
library(gridGraphics)
library(cowplot)

hyperg_combined_info <- read.table("/mfs/gdouglas/projects/ocean_mags/networks/allsamples/clusterbased/metaG_hyperg_cooccur.allsamples.combined.tsv.gz",
                                   header=TRUE, sep="\t", stringsAsFactors = FALSE)
hyperg_combined_info <- hyperg_combined_info[which(! hyperg_combined_info$diff_tax_level %in% c('Species', 'Strain')), ]

hyperg_combined_info$cooccur <- 'No'
hyperg_combined_info$cooccur_ratio <- hyperg_combined_info$cooccur_obs / hyperg_combined_info$cooccur_exp
hyperg_combined_info$cooccur[which(hyperg_combined_info$cooccur_BH < 0.05 & hyperg_combined_info$cooccur_ratio > 1)] <- 'Yes'

hyperg_combined_info$hgt <- 'No'
hyperg_combined_info$hgt[which(hyperg_combined_info$both_gene_count > 0)] <- 'Yes'


freeliving_genomes_file='/mfs/gdouglas/projects/ocean_hgt_zenodo/cooccur/genomes_freeliving_associated.tsv.gz'
lessfiltered_genomes_file='/mfs/gdouglas/projects/ocean_hgt_zenodo/cooccur/genomes_lessfiltered_associated.tsv.gz'
genomes_mainly_freeliving <- read.table(file = freeliving_genomes_file, header=FALSE, stringsAsFactors = FALSE)$V1
genomes_mainly_lessfiltered <- read.table(file = lessfiltered_genomes_file, header=FALSE, stringsAsFactors = FALSE)$V1

taxon_i_grouping <- rep('Other', nrow(hyperg_combined_info))
taxon_j_grouping <- rep('Other', nrow(hyperg_combined_info))
taxon_i_grouping[which(hyperg_combined_info$taxon_i %in% genomes_mainly_freeliving)] <- 'freeliving'
taxon_i_grouping[which(hyperg_combined_info$taxon_i %in% genomes_mainly_lessfiltered)] <- 'lessfiltered'
taxon_j_grouping[which(hyperg_combined_info$taxon_j %in% genomes_mainly_freeliving)] <- 'freeliving'
taxon_j_grouping[which(hyperg_combined_info$taxon_j %in% genomes_mainly_lessfiltered)] <- 'lessfiltered'

hyperg_combined_info$filter_group_match <- 'Mixed'
hyperg_combined_info[which(taxon_i_grouping == 'freeliving' & taxon_j_grouping == 'freeliving'), 'filter_group_match'] <- 'Free-living'
hyperg_combined_info[which(taxon_i_grouping == 'lessfiltered' & taxon_j_grouping == 'lessfiltered'), 'filter_group_match'] <- 'Less-filtered'


table(hyperg_combined_info[which(hyperg_combined_info$filter_group_match == 'Less-filtered'), 'hgt'])
table(hyperg_combined_info[which(hyperg_combined_info$filter_group_match == 'Free-living'), 'hgt'])
table(hyperg_combined_info[which(hyperg_combined_info$filter_group_match == 'Mixed'), 'hgt'])

