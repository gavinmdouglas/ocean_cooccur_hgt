rm(list = ls(all.names = TRUE))

# Create combined table of all BLAST hit-relevant data for use in gene-specific (rather than representative) analyses.

gene_info_w_func <- read.table("/mfs/gdouglas/projects/ocean_mags/functional_annot/gene_info_w_annot.tsv.gz",
                               header=TRUE, sep="\t", stringsAsFactors = FALSE, row.names=1)

COGcategory_and_generep <- gene_info_w_func[, c('COG_category', 'gene_rep')]
COGcategory_and_generep <- COGcategory_and_generep[which(! duplicated(COGcategory_and_generep)), ]

# Then read in actual BLAST output (and taxa identities too).
identities <- c('95', '99')

taxa_breakdown <- read.table('/mfs/gdouglas/projects/ocean_mags/mapfiles/MAG_taxa_breakdown.tsv.gz',
                             header = TRUE, sep = '\t', row.names =1 )

for (identity_cutoff in identities) {

  for (tax_level in c('Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus')) {

    hgt_category <- paste("BLASTn_", tax_level, '_', identity_cutoff, sep = '')

    gene_info_w_func[, hgt_category] <- 'No'

    bed_path <- paste('/mfs/gdouglas/projects/ocean_mags/blast_output/intersections/blast_hit_gene_beds/Hits_',
                      tax_level,
                      '_',
                      identity_cutoff,
                      '_genes_with_B.bed',
                      sep = '')

    bed <- read.table(bed_path, sep = '\t', stringsAsFactors = FALSE, quote = '', comment.char = '')

    # Ignore lines where the gene hit is not present in the mapfile of genes to representatives,
    # as these were excluded from clustering (as they were RNAs)
    unique_gene_hits <- unique(bed$V4[which(bed$V4 %in% rownames(gene_info_w_func))])

    gene_info_w_func[unique_gene_hits, hgt_category] <- 'Yes'

  }

}

clusterbased_hgt <- read.table(file = '/mfs/gdouglas/projects/ocean_mags/clusters/all_best_hits_w_repID.tsv.gz', header = TRUE,
                               sep = '\t', stringsAsFactors = FALSE)

clusterbased_hgt <- clusterbased_hgt[which(clusterbased_hgt$gene1 %in% rownames(gene_info_w_func) & clusterbased_hgt$gene2 %in% rownames(gene_info_w_func)), ]

for (identity_cutoff in identities) {

  if (identity_cutoff == '95') {
    matching_identity_rows_i <- which(clusterbased_hgt$identity >= 95 & clusterbased_hgt$identity < 99)
  } else{
    matching_identity_rows_i <- which(clusterbased_hgt$identity >= 99)
  }

  for (tax_level in c('Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus')) {

    hgt_category <- paste("Clusterbased_", tax_level, '_', identity_cutoff, sep = '')

    gene_info_w_func[, hgt_category] <- 'No'

    matching_taxon_rows_i <- which(clusterbased_hgt$highest_tax_diff == tax_level)

    genes_to_set <- unique(c(clusterbased_hgt$gene1[intersect(matching_identity_rows_i, matching_taxon_rows_i)],
                             clusterbased_hgt$gene2[intersect(matching_identity_rows_i, matching_taxon_rows_i)]))

    gene_info_w_func[genes_to_set, hgt_category] <- 'Yes'

  }

}

analyzed_species <- read.table('/mfs/gdouglas/projects/ocean_mags/species_DTL_analyses/species_to_analyze.txt',
                               stringsAsFactors = FALSE, sep = '', header = FALSE)$V1
raw_data <- list()

for (species in analyzed_species) {

  hgt_bed <- read.table(paste('/mfs/gdouglas/projects/ocean_mags/species_DTL_analyses/gene_beds/', species, '.hgt.bed', sep = ''),
                        stringsAsFactors = FALSE, sep = '\t', header = FALSE)
  hgt_bed_info <- read.table(text = hgt_bed$V4, sep = '|', stringsAsFactors = FALSE)[, -c(2, 4)]
  colnames(hgt_bed_info) <- c('genome', 'gene_family', 'gene')
  hgt_bed_info$species <- species
  hgt_bed_info$hgt <- 'Yes'
  hgt_bed_info$scaffold <- hgt_bed$V1

  nonhgt_bed <- read.table(paste('/mfs/gdouglas/projects/ocean_mags/species_DTL_analyses/gene_beds/', species, '.nonhgt.bed', sep = ''),
                           stringsAsFactors = FALSE, sep = '\t', header = FALSE)
  nonhgt_bed_info <- read.table(text = nonhgt_bed$V4, sep = '|', stringsAsFactors = FALSE)[, -c(2, 4)]
  colnames(nonhgt_bed_info) <- c('genome', 'gene_family', 'gene')
  nonhgt_bed_info$species <- species
  nonhgt_bed_info$hgt <- 'No'
  nonhgt_bed_info$scaffold <- nonhgt_bed$V1

  raw_data[[species]] <- rbind(hgt_bed_info, nonhgt_bed_info)
}

ranger_data <- do.call(rbind, raw_data)
rownames(ranger_data) <- ranger_data$gene

ranger_data <- ranger_data[which(ranger_data$gene %in% rownames(gene_info_w_func)), ]

gene_info_w_func$ranger_species <- NA
gene_info_w_func$ranger_hgt <- NA

gene_info_w_func[ranger_data$gene, 'ranger_species'] <- ranger_data$species
gene_info_w_func[ranger_data$gene, 'ranger_hgt'] <- ranger_data$hgt

genes_to_panaroo <- read.table(file = "/mfs/gdouglas/projects/ocean_mags/species_DTL_analyses/genes_to_panaroo_families.tsv.gz",
                               sep = "\t", stringsAsFactors = FALSE, header = TRUE )
genes_to_panaroo <- genes_to_panaroo[-grep("refound", genes_to_panaroo$gene_id), ]
rownames(genes_to_panaroo) <- genes_to_panaroo$gene_id

intersecting_genes <- intersect(genes_to_panaroo$gene_id, rownames(gene_info_w_func))

gene_info_w_func[intersecting_genes, "panaroo_sp_gf"] <- genes_to_panaroo[intersecting_genes, "panaroo_sp_and_gene_family"]

orig_col <- colnames(gene_info_w_func)

gene_info_w_func$gene <- rownames(gene_info_w_func)

gene_info_w_func <- gene_info_w_func[, c("gene", orig_col)]

write.table(x = gene_info_w_func,
            file = '/mfs/gdouglas/projects/ocean_mags/summary_files/gene_info_and_hgt_redone2.tsv',
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = '\t')
