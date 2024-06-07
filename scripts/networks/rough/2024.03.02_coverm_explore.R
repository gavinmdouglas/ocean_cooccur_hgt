rm(list = ls(all.names = TRUE))

bioprojects <- c('PRJEB1787', 'PRJEB6608', 'PRJEB97740')

raw_in <- list()

for (bioproject in bioprojects) {
 
  all_files <- list.files(paste0("/Users/gavin/Drive/mcgill/water_mags_analyses/coverm_output/", bioproject), 
                          pattern = "cov.gz", full.names = TRUE)
  
  cat("Bioproject:", bioproject, "\n")
  cat("Number of samples:", length(all_files), "\n")
  cat("\n\n")
  
  for (file in all_files) {
    cov_in <- read.table(file, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
    
    colnames(cov_in) <- c('genome', 'mean', 'rpkm', 'trimmed_mean', 'relabun', 'covered_bases', 'length')

    cov_in$breadth <- cov_in$covered_bases / cov_in$length
    
    orig_col <- colnames(cov_in)
    
    mgs_id <- gsub('\\..*$', '', basename(file))
    cov_in$mgs_sample <- mgs_id

    cov_in$bioproject <- bioproject
    
    cov_in <- cov_in[, c('bioproject', 'mgs_sample', orig_col)]

    raw_in[[paste(bioproject, mgs_id)]] <- cov_in
  }
   
}

combined <- do.call(rbind, raw_in)
rownames(combined) <- NULL

# Sanity check that all mgs samples have the same number of rows.
# (Since all genomes mapped to should be represented, even if no reads mapepd).
table(table(combined$mgs_sample))

# Add in taxonomic info.
taxa_map <- read.table('/Users/gavin/Drive/mcgill/water_mags_analyses/MAG_taxa_breakdown.tsv.gz',
                       header = TRUE, sep = '\t', stringsAsFactors = FALSE, row.names = 2)

species_tallies <- table(taxa_map$Species)
species_breakdown <- data.frame(count = as.integer(species_tallies))
rownames(species_breakdown) <- names(species_tallies)

combined$species <- taxa_map[combined$genome, 'Species']

# Call taxa as present (1) if > 10% of the genome is covered, otherwise set 0.
combined$present0.1 <- ifelse(combined$breadth > 0.1, 1, 0)
combined$present0.3 <- ifelse(combined$breadth > 0.3, 1, 0)

test_sp <- "d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Pelagibacterales; f__Pelagibacteraceae; g__Pelagibacter; s__"

combined_subset <- combined[which(combined$species == test_sp), ]

combined_subset_num_present0.1 <- aggregate(present0.1 ~ mgs_sample, data = combined_subset, FUN = sum)
combined_subset_num_present0.3 <- aggregate(present0.3 ~ mgs_sample, data = combined_subset, FUN = sum)

combined$present0.75 <- ifelse(combined$breadth > 0.75, 1, 0)
num_genomes_per_sample0.75 <- aggregate(present0.75 ~ mgs_sample, data = combined, FUN = sum)

combined$present0.5 <- ifelse(combined$relabun > 0.5, 1, 0)
num_genomes_per_sample0.5 <- aggregate(present0.5 ~ mgs_sample, data = combined, FUN = sum)


# In this case, ignore cases where no genomes for the species are present.
combined_subset_num_present0.1 <- combined_subset_num_present0.1[which(combined_subset_num_present0.1$present0.1 > 0), ]
combined_subset_num_present0.3 <- combined_subset_num_present0.3[which(combined_subset_num_present0.3$present0.3 > 0), ]

# Get number of genomes per sample.
num_genomes_per_sample0.1 <- aggregate(present0.1 ~ mgs_sample, data = combined, FUN = sum)
num_genomes_per_sample0.3 <- aggregate(present0.3 ~ mgs_sample, data = combined, FUN = sum)


# Dig into whether the number of genomes per species is negatively associated with the max breadth
# (which you might expect if there is too much mapping competition).
combined_above0.1 <- combined[which(combined$breadth > 0.1), ]
combined_above0.1_max_per_species_sample <- aggregate(breadth ~ mgs_sample + species, data = combined_above0.1, FUN = max)

combined_above0.1_max_per_averaged <- aggregate(breadth ~ species, data = combined_above0.1_max_per_species_sample, FUN = mean)

# Add numbers of genomes in each species.
combined_above0.1_max_per_averaged$species_count <- species_breakdown[combined_above0.1_max_per_averaged$species, 'count']
plot(combined_above0.1_max_per_averaged$species_count, combined_above0.1_max_per_averaged$breadth)
cor.test(combined_above0.1_max_per_averaged$species_count, combined_above0.1_max_per_averaged$breadth)

