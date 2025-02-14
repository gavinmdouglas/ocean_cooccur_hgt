rm(list = ls(all.names = TRUE))

library(taxize)

# packageVersion('taxize')
# [1] ‘0.9.102’

taxids <- read.table('/mfs/gdouglas/projects/ocean_mags/progenomes_analyses/genome_prep/all_taxids.txt', stringsAsFactors = FALSE, header=FALSE)$V1

taxonomy_raw <- lapply(taxids, function(x) classification(x, db = "ncbi", return_id = TRUE))

saveRDS(object = taxonomy_raw, file = "/mfs/gdouglas/projects/ocean_mags/progenomes_analyses/genome_prep/taxid_taxize_download.rds")
