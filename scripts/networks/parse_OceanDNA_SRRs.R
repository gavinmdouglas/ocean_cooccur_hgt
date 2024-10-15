rm(list = ls(all.names = TRUE))

# Parse OceanDNA SRRs to those that we need to download.
# Also check that this matches the set that Nico downloaded previously.
# Then output file with all SRRs to download, split into separate BioProjects.

library("readxl")

TableS1 <- read_excel("/mfs/gdouglas/projects/ocean_mags/metadata/OceanDNA_supp_metadata/Supp_File_S1.xlsx", skip=2)

Nico_subset <- read_excel("/mfs/gdouglas/projects/ocean_mags/metadata/OceanDNA_supp_metadata/Samples_to_keep.xlsx", skip=1)

orig_bioprojects <- unique(TableS1$bioproject)
subset_bioprojects <- unique(Nico_subset$bioproject)
# Sanity checked that all removed bioprojects made sense.
# Note sure about PRJNA438384, but filters were set to very large particles, so seems reasonable to exclude.

for (bioproject in subset_bioprojects) {
     
  tab_subset <- Nico_subset[which(Nico_subset$bioproject == bioproject), ]
  
  raw_run_ids <- tab_subset$sra_run
  run_ids <- character()
  for (i in 1:length(raw_run_ids)) {
    run_ids <- c(run_ids, strsplit(raw_run_ids[i], ",")[[1]])
  }
  
  outfile <- paste0("/mfs/gdouglas/projects/ocean_mags/additional/id_by_bioproject/",
                   bioproject, ".txt")
  
  write.table(x = run_ids, file = outfile, quote = FALSE, col.names = FALSE, row.names = FALSE)
    
}

# Also get simple table in TSV format (for subset of samples of interest).
write.table(x = Nico_subset,
            file = "/mfs/gdouglas/projects/ocean_mags/additional/OceanDNA_supp_metadata/subset_tab.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
