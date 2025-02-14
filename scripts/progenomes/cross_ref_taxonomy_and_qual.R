rm(list = ls(all.names = TRUE))

# Cross-reference taxonomy breakdown with list of high quality genomes.

# Can only keep genomes (as "high quality") that are defined at least at genus level.

taxa <- read.table('/mfs/gdouglas/projects/ocean_mags/progenomes_analyses/genome_prep/progenomes_ncbi_taxonomy.tsv.gz',
                   header = TRUE, stringsAsFactors = FALSE, sep = '\t')

highqual_ids <- read.table('/mfs/gdouglas/projects/ocean_mags/progenomes_analyses/genome_prep/clearly_non_mags_genome_ids.txt',
                           header = FALSE, stringsAsFactors = FALSE, sep = '\t')$V1

highqual_taxids <- gsub('_.*$', '', highqual_ids)
highqual_taxids <- gsub('^g', '', highqual_taxids)

# Get taxa table with duplicate by biosample.
rownames(taxa) <- taxa$taxid
taxa_all_highqual <- taxa[highqual_taxids, ]

taxa_all_highqual$taxid_biosample <- highqual_ids

taxa_all_highqual_defined_genus <- taxa_all_highqual[grep(' NA$',taxa_all_highqual$Genus, invert = TRUE), ]

write.table(file = '/mfs/gdouglas/projects/ocean_mags/progenomes_analyses/genome_prep/clearly_non_mags_nonambig.genus_genome_ids.txt',
            x = taxa_all_highqual_defined_genus$taxid_biosample,
            col.names = FALSE, row.names = FALSE, quote = FALSE)
