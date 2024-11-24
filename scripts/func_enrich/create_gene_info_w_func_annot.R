rm(list = ls(all.names = TRUE))

# Create combined table of all gene info with functional annotations.
# Subset to genes on scaffold >= 5000 bp only.

gene_info <- read.table(file = '/mfs/gdouglas/projects/ocean_mags/Sunagawa_dataset/gene_info.tsv.gz',
                        sep = '\t', stringsAsFactors = FALSE, header = TRUE)
rownames(gene_info) <- gene_info$gene

# Read in COG categories
COG_annot <- read.table('/mfs/gdouglas/projects/ocean_mags/Sunagawa_dataset/genomes-cog-info.tsv.gz',
                        header = TRUE, sep = '\t', stringsAsFactors = FALSE)
rownames(COG_annot) <- COG_annot$query_name

# Set all genes missing COG annotation as '-'.
gene_info$COG_category <- '-'
# Set all others to be what was specified.
COG_annot_genes <- intersect(rownames(gene_info), rownames(COG_annot))
gene_info[COG_annot_genes, 'COG_category'] <- COG_annot[COG_annot_genes, 'COG_category']

# Then specify endonuclease/recombinase annotation as variable.
proMGE_annot <- read.table('/mfs/gdouglas/projects/ocean_mags/functional_annot/recombinanse_and_endonucleases/Sunagawa_proMGE_v1_hmmersearch_tab.edit.tsv.gz',
                           header = TRUE, stringsAsFactors = FALSE)
proMGE_annot$target_name <- gsub('_1$', '', proMGE_annot$target_name)

gene_info$proMGE <- 'No'
gene_info$proMGE[which(gene_info$gene %in% proMGE_annot$target_name)] <- 'Yes'


# Add in geNomad scores for chromosome/plasmid/virus.
genomad_virus_scaffolds <- read.table('/mfs/gdouglas/projects/ocean_mags/functional_annot/genomad_out/genomad_out_MAGS/noDuplicates_virus_contigs.txt',
                                      header = FALSE, stringsAsFactors = FALSE)$V1
genomad_plasmid_scaffolds <- read.table('/mfs/gdouglas/projects/ocean_mags/functional_annot/genomad_out/genomad_out_MAGS/noDuplicates_plasmids_contigs.txt',
                                        header = FALSE, stringsAsFactors = FALSE)$V1

gene_info$genomad <- 'Other'
gene_info$genomad[which(gene_info$scaffold %in% genomad_virus_scaffolds)] <- 'Virus'
gene_info$genomad[which(gene_info$scaffold %in% genomad_plasmid_scaffolds)] <- 'Plasmid'

# Mark scaffolds containing proviruses.
provirus_breakdown <- read.table('/mfs/gdouglas/projects/ocean_mags/functional_annot/genomad_out/provirus_breakdown.tsv.gz',
                                 header = FALSE, stringsAsFactors = FALSE, sep = '\t')
gene_info$scaffold_contains_provirus <- 'No'
gene_info[which(gene_info$scaffold %in% provirus_breakdown$V1), "scaffold_contains_provirus"] <- 'Yes'

# Make sure only genes on scaffolds that are at least 5000 bp are included here.
scaffolds_to_keep <- read.table("/mfs/gdouglas/projects/ocean_mags/Sunagawa_dataset/scaffolds_5000bp.txt", stringsAsFactors = FALSE, header=FALSE)$V1

gene_info <- gene_info[which(gene_info$scaffold %in% scaffolds_to_keep), ]

write.table(x = gene_info,
            file = '/mfs/gdouglas/projects/ocean_mags/functional_annot/gene_info_w_annot.tsv',
            sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)
