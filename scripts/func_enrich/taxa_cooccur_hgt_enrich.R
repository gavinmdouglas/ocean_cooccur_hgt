rm(list = ls(all.names = TRUE))

# Identify taxa enriched as being especially more (or less) often co-occurring
# when they undergo HGT.

library(ggplot2)

taxa_cooccur_hgt <- read.table("/mfs/gdouglas/projects/ocean_mags/coverm/network_working_allsamples/prepped_taxa_cooccur_vs_nonoccur_hgt.tsv.gz",
                               header = TRUE, sep = "\t", stringsAsFactors = FALSE)

taxa_cooccur_hgt$total_hgt <- taxa_cooccur_hgt$Cooccur_HGT + taxa_cooccur_hgt$No_Cooccur_HGT

taxa_cooccur_hgt_atleast20 <- taxa_cooccur_hgt[which(taxa_cooccur_hgt$total_hgt >= 20), ]

taxa_cooccur_hgt$Tax_Level <- factor(taxa_cooccur_hgt$Tax_Level,
                                     levels = c('Genus', 'Family', 'Order', 'Class', 'Phylum', 'Domain'))

fisher_out <- lapply(1:nrow(taxa_cooccur_hgt_atleast20),
                                        function(row_i) {
                                          
                                          num_coccur_hgt <- taxa_cooccur_hgt_atleast20[row_i, 'Cooccur_HGT']
                                          num_coccur_nohgt <- taxa_cooccur_hgt_atleast20[row_i, 'Cooccur_No_HGT']
                                          num_nococcur_hgt <- taxa_cooccur_hgt_atleast20[row_i, 'No_Cooccur_HGT']
                                          num_nococcur_nohgt <- taxa_cooccur_hgt_atleast20[row_i, 'No_Cooccur_No_HGT']
                                          
                                          OR <- ((num_coccur_hgt + 1) / (num_coccur_nohgt + 1)) / ((num_nococcur_hgt + 1) / (num_nococcur_nohgt + 1))
                                          
                                          P <- fisher.test(matrix(c(num_coccur_hgt, num_coccur_nohgt, num_nococcur_hgt, num_nococcur_nohgt), nrow = 2))$p.value
                                          
                                          return(list(OR=OR, P=P))
                                          
                                        })

taxa_cooccur_hgt_atleast20$OR <- sapply(fisher_out, function(x) { x$OR })
taxa_cooccur_hgt_atleast20$P <- sapply(fisher_out, function(x) { x$P })

# Keep only one domain taxon, as they are 100% redundant.
taxa_cooccur_hgt_atleast20 <- taxa_cooccur_hgt_atleast20[-which(taxa_cooccur_hgt_atleast20$Full_Taxonomy == "d__Archaea"), ]

taxa_cooccur_hgt_atleast20$BH <- p.adjust(taxa_cooccur_hgt_atleast20$P, 'BH')

taxa_cooccur_hgt_atleast20 <- taxa_cooccur_hgt_atleast20[order(taxa_cooccur_hgt_atleast20$OR, decreasing = TRUE), ]
taxa_cooccur_hgt_atleast20$Full_Taxonomy <- factor(taxa_cooccur_hgt_atleast20$Full_Taxonomy, levels = taxa_cooccur_hgt_atleast20$Full_Taxonomy)

y_breaks <- taxa_cooccur_hgt_atleast20$Full_Taxonomy[seq(1, length(taxa_cooccur_hgt_atleast20$Full_Taxonomy), by = 3)]

ggplot(data = taxa_cooccur_hgt_atleast20,
       aes(x = log2(OR), y = Full_Taxonomy, fill = Tax_Level)) +
  geom_col() +
  scale_y_discrete(breaks = y_breaks) +
  theme_bw() +
  ylab('Taxon') +
  xlab(expression(log[2]("Odd's ratio"))) +
  labs(fill="Tax. level") +
  scale_fill_manual(values=c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02"))
