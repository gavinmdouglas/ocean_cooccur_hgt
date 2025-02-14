rm(list = ls(all.names = TRUE))

info <- read.table('/mfs/gdouglas/projects/ocean_mags/progenomes_analyses/genome_prep/metagenomic_and_status_info.tsv',
                   header=TRUE, stringsAsFactors = FALSE, sep = '\t', row.names=1, quote="", comment.char = "")

convert_str_percent_to_num <- function(x) {
  if (is.na(x)) {
    return(NA)
  } else {
    x <- gsub('%', '', x)
    x <- as.numeric(gsub(',', '.', x))
    return(x)
  }

}

info$contam_clean <- NA
info$complete_clean <- NA

for (i in 1:nrow(info)) {
  contam1 <- convert_str_percent_to_num(info$contamination_estimated[i])
  contam2 <- convert_str_percent_to_num(info$contamination[i])
  contam3 <- convert_str_percent_to_num(info$contamination.score[i])
  contam_values <- unique(c(contam1, contam2, contam3))

  if (length(which(is.na(contam_values))) > 0) {
    contam_values <- contam_values[-which(is.na(contam_values))]
  }

  if (length(contam_values) == 1) {
    info$contam_clean[i] <- contam_values[1]
  } else if (length(contam_values) > 1) {
    stop('More than 1 contamination value for row ', i)
  }

  complete1 <- convert_str_percent_to_num(info$completeness_estimated[i])
  complete2 <- convert_str_percent_to_num(info$completeness[i])
  complete3 <- convert_str_percent_to_num(info$completeness.score[i])
  complete_values <- unique(c(complete1, complete2, complete3))

  if (length(which(is.na(complete_values))) > 0) {
    complete_values <- complete_values[-which(is.na(complete_values))]
  }

  if (length(complete_values) == 1) {
    info$complete_clean[i] <- complete_values[1]
  } else if (length(complete_values) > 1) {
    stop('More than 1 completeness value for row ', i)
  }
}

tallies <- read.table('/mfs/gdouglas/projects/ocean_mags/progenomes_analyses/genome_prep/contig_or_chrom_per_genome.txt',
                      header=TRUE, stringsAsFactors = FALSE)

tallies$biosample <- gsub("^.*\\.", "", tallies$genome)

tallies$sample_type <- info[tallies$biosample, 'sample.type']
tallies$isolation_source <- info[tallies$biosample, 'isolation.source']
tallies$isolate <- info[tallies$biosample, 'isolate']
tallies$metagenome_source <- info[tallies$biosample, 'metagenome.source']
tallies$metagenomic <- info[tallies$biosample, 'metagenomic']
tallies$complete_clean <- info[tallies$biosample, 'complete_clean']
tallies$contam_clean <- info[tallies$biosample, 'contam_clean']

tallies_no_mags <- tallies[which(is.na(tallies$metagenomic)), ]
tallies_no_mags <- tallies_no_mags[which(is.na(tallies_no_mags$metagenome_source)), ]

tallies_no_mags <- tallies_no_mags[which(tallies_no_mags$MAG_in_name == 'No'), ]
tallies_no_mags <- tallies_no_mags[grep('metagenome', tallies_no_mags$isolation_source, invert = TRUE), ]

# Get this set written out.
non_mag_genomes <- paste('g', tallies_no_mags$genome, sep = '')
non_mag_genomes <- gsub('\\.', '_', non_mag_genomes)

# Also subset to those with at least one contig above length cut-off.
contig_lengths <- read.table('/mfs/gdouglas/projects/ocean_mags/progenomes_analyses/genome_prep/contig_length_breakdown.tsv',
                             header=TRUE, sep = '\t', stringsAsFactors = FALSE)

contig_lengths_long <- contig_lengths[which(contig_lengths$type == 'long'), ]

non_mag_genomes_nonshort <- non_mag_genomes[which(non_mag_genomes %in% contig_lengths_long$genome)]

write.table(x = non_mag_genomes_nonshort, file = '/mfs/gdouglas/projects/ocean_mags/progenomes_analyses/genome_prep/clearly_non_mags_genome_ids.txt',
            quote=FALSE, row.names = FALSE, col.names = FALSE)



# Also get tallies broken down by plasmid vs other (assumed to be chromosomal) contig.
plasmid_vs_other <- read.table('/mfs/gdouglas/projects/ocean_mags/progenomes_analyses/genome_prep/contig_tally_plasmid_vs_other.tsv',
                               header=TRUE, stringsAsFactors = FALSE, row.names = 1)


# # Get genomes with single chrom in FASTA (besides those marked as plasmids).
# single_chrom_genomes <- read.table('/mfs/gdouglas/projects/ocean_mags/progenomes_analyses/genome_prep/single_chrom_genomes.txt',
#                                    header=FALSE, stringsAsFactors = FALSE)$V1
#
# tallies_single_chrom <- tallies[tallies$genome %in% single_chrom_genomes, ]
#
# tallies_single_chrom_no_mags <- tallies_single_chrom[which(is.na(tallies_single_chrom$metagenomic)), ]
# tallies_single_chrom_no_mags <- tallies_single_chrom_no_mags[which(is.na(tallies_single_chrom_no_mags$metagenome_source)), ]
#
# tallies_single_chrom_no_mags <- tallies_single_chrom_no_mags[which(tallies_single_chrom_no_mags$MAG_in_name == 'No'), ]
# tallies_single_chrom_no_mags <- tallies_single_chrom_no_mags[grep('metagenome', tallies_single_chrom_no_mags$isolation_source, invert = TRUE), ]
#
# write.table(x = tallies_single_chrom_no_mags$genome, file = '/mfs/gdouglas/projects/ocean_mags/progenomes_analyses/genome_prep/single_chrom_non_mags_genome_ids.txt',
#             quote=FALSE, row.names = FALSE, col.names = FALSE)
