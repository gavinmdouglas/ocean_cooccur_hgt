rm(list = ls(all.names = TRUE))

# Quick double-check that focal genome filtering was based on > 75% completeness and < 5% contamination.
genome_meta <- read.table('~/projects/ocean_mags/Sunagawa_dataset/genomes-summary.csv.gz',
                          header=TRUE, sep = ',', stringsAsFactors = FALSE, row.names=1)

retained_genomes <- read.table('~/projects/ocean_mags/mapfiles/MAGs_to_analyze.txt.gz', stringsAsFactors = FALSE, header=FALSE)$V1

setdiff(retained_genomes, rownames(genome_meta))

retained_meta <- genome_meta[retained_genomes, ]

summary(retained_meta$Mean.Completeness)
summary(retained_meta$Mean.Contamination)

# Number total
nrow(genome_meta)
# [1] 34815

# Number retained
length(retained_genomes)
# [1] 15339


# Get counts of genomes by sequencing method/source:
mag_ids <- grep("_METAG_", rownames(retained_meta), value=TRUE)
sag_ids <- grep("_SAGS_", rownames(retained_meta), value=TRUE)
ref_ids <- grep("_REFG_", rownames(retained_meta), value=TRUE)

length(mag_ids)
length(sag_ids)
length(ref_ids)

# Sanity check that all genomes only match 1 category, and that no genomes missing.
mag_sag_overlap <- intersect(mag_ids, sag_ids)
mag_ref_overlap <- intersect(mag_ids, ref_ids)
sag_ref_overlap <- intersect(sag_ids, ref_ids)

if(length(mag_sag_overlap) > 0 | length(mag_ref_overlap) > 0 | length(sag_ref_overlap) > 0) {
  warning("Overlaps detected between categories!")
}

all_categorized <- c(mag_ids, sag_ids, ref_ids)
uncategorized <- setdiff(rownames(retained_meta), all_categorized)
if(length(uncategorized) > 0) {
  cat("\nUncategorized IDs (first 10):\n")
  print(head(uncategorized, 10))
}
