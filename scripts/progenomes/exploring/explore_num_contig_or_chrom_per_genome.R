rm(list = ls(all.names = TRUE))

library(ggplot2)
library(reshape2)
library(ComplexHeatmap)
library(gridGraphics)
library(cowplot)

info <- read.table('/mfs/gdouglas/projects/ocean_mags/progenomes_analyses/genome_prep/metagenomic_and_status_info.tsv',
                   header=TRUE, stringsAsFactors = FALSE, sep = '\t', row.names=1, quote="", comment.char = "")

# All present except for 6 wonky IDs!

# Two suppressed:
# 5955 SAMN09639835 Present <NA> suppressed
# 6689 SAMN16057880 Present <NA> suppressed

tallies <- read.table('/mfs/gdouglas/projects/ocean_mags/progenomes_analyses/genome_prep/contig_or_chrom_per_genome.txt',
                      header=TRUE, stringsAsFactors = FALSE)

tallies$biosample <- gsub("^.*\\.", "", tallies$genome)

tallies$sample_type <- info[tallies$biosample, 'sample.type']
tallies$isolation_source <- info[tallies$biosample, 'isolation.source']
tallies$isolate <- info[tallies$biosample, 'isolate']
tallies$metagenome_source <- info[tallies$biosample, 'metagenome.source']
tallies$metagenomic <- info[tallies$biosample, 'metagenomic']


# Get genomes with single chrom in FASTA (besides those marked as plasmids).
single_chrom_genomes <- read.table('/mfs/gdouglas/projects/ocean_mags/progenomes_analyses/genome_prep/single_chrom_genomes.txt',
                                   header=FALSE, stringsAsFactors = FALSE)$V1

tallies_single_chrom <- tallies[tallies$genome %in% single_chrom_genomes, ]

tallies_single_chrom_no_mags <- tallies_single_chrom[which(is.na(tallies_single_chrom$metagenomic)), ]
tallies_single_chrom_no_mags <- tallies_single_chrom_no_mags[which(is.na(tallies_single_chrom_no_mags$metagenome_source)), ]

tallies_single_chrom_no_mags <- tallies_single_chrom_no_mags[which(tallies_single_chrom_no_mags$MAG_in_name == 'No'), ]
tallies_single_chrom_no_mags <- tallies_single_chrom_no_mags[grep('metagenome', tallies_single_chrom_no_mags$isolation_source, invert = TRUE), ]

tallies_single_chrom_no_mags$

tallies$mgs_sample_type <- FALSE
tallies[grep('metagenom', tallies$sample_type), 'mgs_sample_type'] <- TRUE
tallies[grep('Metagenomic', tallies$sample_type), 'mgs_sample_type'] <- TRUE

tallies_nonmgs <- tallies[which(! tallies$mgs_sample_type), ]
tallies_nonmgs <- tallies_nonmgs[which(tallies_nonmgs$contig_count <= 10), ]


tallies$any_MAG_in_name <- 'No'
tallies[which(tallies$MAG_in_name == 'Yes' | tallies$mag_in_name == 'Yes'), 'any_MAG_in_name'] <- 'Yes'

ggplot(data = tallies, aes(x = any_MAG_in_name, y=contig_count )) +
  geom_boxplot() +
  theme_bw()

mag_tallies <- tallies[which(tallies$any_MAG_in_name == 'Yes'), ]
nonmag_tallies <- tallies[which(tallies$any_MAG_in_name == 'No'), ]


tallies_filt <- tallies[w]