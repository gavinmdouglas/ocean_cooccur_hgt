rm(list = ls(all.names = TRUE))

# Quick check that focal genome filtering was based on > 75% completeness and < 5% contamination.

genome_meta <- read.table('~/projects/ocean_mags/Sunagawa_dataset/genomes-summary.csv.gz',
                          header=TRUE, sep = ',', stringsAsFactors = FALSE, row.names=1)

retained_genomes <- read.table('~/projects/ocean_mags/mapfiles/MAGs_to_analyze.txt.gz', stringsAsFactors = FALSE, header=FALSE)$V1

setdiff(retained_genomes, rownames(genome_meta))

retained_meta <- genome_meta[retained_genomes, ]

hist(retained_meta$Mean.Completeness)

hist(retained_meta$Mean.Contamination)

