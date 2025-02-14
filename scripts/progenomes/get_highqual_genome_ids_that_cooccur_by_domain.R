rm(list = ls(all.names = TRUE))

tmp <- read.table('/mfs/gdouglas/projects/ocean_mags/progenomes_analyses/networks/allsamples/clusterbased/metaG_hyperg_cooccur.allsamples.combined.tsv', sep = '\t', stringsAsFactors = FALSE, header = TRUE)
length(which(tmp$cooccur_BH < 0.05 & tmp$cooccur_ratio > 1))

length(which(tmp$cooccur_BH < 0.05 & tmp$cooccur_ratio > 1 & tmp$both_gene_count > 0))
length(which(tmp$cooccur_BH < 0.05 & tmp$cooccur_ratio <= 1 & tmp$both_gene_count > 0))
length(which(tmp$cooccur_BH > 0.05 & tmp$both_gene_count > 0))
length(which(tmp$cooccur_BH >= 0.05 & tmp$both_gene_count > 0))
length(which(tmp$cooccur_BH >= 0.05))


all_highqual <- unique(c(tmp$taxon_i, tmp$taxon_j))

all_highqual_taxid <- sapply(all_highqual,
                             function(x) {
    x <- gsub('^g', '', x)
    strsplit(x = x, split = '_')[[1]][1]
})


genome_taxonomy <- read.table('/mfs/gdouglas/projects/ocean_mags/progenomes_analyses/genome_prep/progenomes_ncbi_taxonomy.tsv.gz',
                              header=TRUE, sep = '\t', stringsAsFactors = FALSE, row.names=1)

all_highqual_domain <- genome_taxonomy[all_highqual_taxid, 'Domain']

bacteria_highqual_ids <- all_highqual[which(all_highqual_domain == 'Bacteria')]
archaea_highqual_ids <- all_highqual[which(all_highqual_domain == 'Archaea')]

write.table(x = bacteria_highqual_ids, file = '/mfs/gdouglas/projects/ocean_mags/progenomes_analyses/genome_prep/highqual_genome_fastas/bacteria_ids.txt',
            quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)

write.table(x = archaea_highqual_ids, file = '/mfs/gdouglas/projects/ocean_mags/progenomes_analyses/genome_prep/highqual_genome_fastas/archaea_ids.txt',
            quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)


# orig_mags <- read.table('/mfs/gdouglas/projects/ocean_mags/networks/allsamples/clusterbased/metaG_hyperg_cooccur.allsamples.combined.tsv.gz',
#                         sep = '\t', stringsAsFactors = FALSE, header = TRUE)
#
# length(which(orig_mags$cooccur_BH < 0.05 & orig_mags$cooccur_ratio > 1))
# length(which(orig_mags$cooccur_BH < 0.05 & orig_mags$cooccur_ratio > 1 & orig_mags$both_gene_count > 0))
# length(which(orig_mags$cooccur_BH < 0.05 & orig_mags$cooccur_ratio <= 1 & orig_mags$both_gene_count > 0))
# length(which(orig_mags$cooccur_BH > 0.05 & orig_mags$both_gene_count > 0))
# length(which(orig_mags$cooccur_BH >= 0.05 & orig_mags$both_gene_count > 0))
# length(which(orig_mags$cooccur_BH >= 0.05))
