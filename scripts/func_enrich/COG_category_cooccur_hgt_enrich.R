rm(list = ls(all.names = TRUE))

# Identify COG categories HGT events by whether they are more or less likely to
# be between co-occurring genomes.
library(ComplexHeatmap)


func_cooccur_hgt <- read.table("/mfs/gdouglas/projects/ocean_mags/coverm/network_working_allsamples/prepped_COG_cooccur_vs_nonoccur_hgt.tsv.gz",
                               header = TRUE, sep = "\t", stringsAsFactors = FALSE)

to_ignore <- c("A", "B", "Y", "Z")

func_cooccur_hgt <- func_cooccur_hgt[which(! func_cooccur_hgt$COG_Category %in% to_ignore), ]

func_cooccur_hgt$total <- func_cooccur_hgt$Cooccur_HGT_Gene_Count + func_cooccur_hgt$Noncooccur_HGT_Gene_Count

hist(func_cooccur_hgt$total)

func_cooccur_hgt_atleast50 <- func_cooccur_hgt[which(func_cooccur_hgt$total >= 50), ]

func_cooccur_hgt_atleast50$Identity <- as.character(func_cooccur_hgt_atleast50$Identity)

raw_out <- list()

for (identity_cutoff in unique(func_cooccur_hgt_atleast50$Identity)) {
  identity_subset <- func_cooccur_hgt_atleast50[which(func_cooccur_hgt_atleast50$Identity == identity_cutoff), , drop = FALSE]

  for (tax in unique(identity_subset$Tax_Level)) {
    tax_identity_subset <- identity_subset[which(identity_subset$Tax_Level == tax), ]

    subset_cooccur_total <- sum(tax_identity_subset$Cooccur_HGT_Gene_Count)
    subset_nococcur_total <- sum(tax_identity_subset$Noncooccur_HGT_Gene_Count)
  
    fisher_out <- lapply(1:nrow(tax_identity_subset),
                         function(row_i) {
                           
                           num_coccur_hgt <- tax_identity_subset[row_i, 'Cooccur_HGT_Gene_Count']
                           num_nococcur_hgt <- tax_identity_subset[row_i, 'Noncooccur_HGT_Gene_Count']
                           other_coccur_hgt <- subset_cooccur_total - tax_identity_subset[row_i, 'Cooccur_HGT_Gene_Count']
                           other_nococcur_hgt <- subset_nococcur_total - tax_identity_subset[row_i, 'Noncooccur_HGT_Gene_Count']
                           
                           OR <- ((num_coccur_hgt + 1) / (num_nococcur_hgt + 1)) / ((other_coccur_hgt + 1) / (other_nococcur_hgt + 1))
                           
                           P <- fisher.test(matrix(c(num_coccur_hgt, num_nococcur_hgt, other_coccur_hgt, other_nococcur_hgt), nrow = 2))$p.value
                           
                           return(list(num_coccur_hgt=num_coccur_hgt,
                                       num_nococcur_hgt=num_nococcur_hgt,
                                       other_coccur_hgt=other_coccur_hgt,
                                       other_nococcur_hgt=other_nococcur_hgt,
                                       OR=OR,
                                       P=P))
                           
                         })
    
    raw_out[[paste(identity_cutoff, tax)]] <- data.frame(identity = identity_cutoff,
                                                 tax_level = tax,
                                                 category = tax_identity_subset$COG_Category,
                                                 num_coccur_hgt=sapply(fisher_out, function(x) { return(x$num_coccur_hgt)}),
                                                 num_nococcur_hgt=sapply(fisher_out, function(x) { return(x$num_nococcur_hgt)}),
                                                 other_coccur_hgt=sapply(fisher_out, function(x) { return(x$other_coccur_hgt)}),
                                                 other_nococcur_hgt=sapply(fisher_out, function(x) { return(x$other_nococcur_hgt)}),
                                                 OR=sapply(fisher_out, function(x) { return(x$OR)}),
                                                 P=sapply(fisher_out, function(x) { return(x$P)}))
        
  }
}

combined <- do.call(rbind, raw_out)

combined$BH <- p.adjust(combined$P, 'BH')

combined$log2OR <- log2(combined$OR)

combined$tax_level <- factor(combined$tax_level,
                             levels = rev(c('Genus', 'Family', 'Order', 'Class', 'Phylum', 'Domain')))

combined$category <- factor(combined$category)

combined_95 <- combined[which(combined$identity == '95'), ]

combined_95_OR <- reshape2::dcast(data = combined_95, formula = tax_level ~ category, value.var = "log2OR", fill = NA, drop = FALSE)
combined_95_BH <- reshape2::dcast(data = combined_95, formula = tax_level ~ category, value.var = "BH", fill = 1, drop = FALSE)

rownames(combined_95_OR) <- combined_95_OR$tax_level
combined_95_OR <- combined_95_OR[, -1]

rownames(combined_95_BH) <- combined_95_BH$tax_level
combined_95_BH <- combined_95_BH[, -1]

combined_95_OR_char <- apply(combined_95_OR, 2, function(x) { gsub(" *", "", sprintf("%.2f", x)) })
combined_95_OR_char <- apply(combined_95_OR_char, 2, function(x) { gsub(" *", "", x) })
combined_95_OR_char <- apply(combined_95_OR_char, 2, function(x) { gsub("NA", "", x) })
combined_95_OR_char[is.na(combined_95_OR_char)] <- ''

combined_95_OR[combined_95_BH >= 0.05] <- NA


ComplexHeatmap::Heatmap(matrix = as.matrix(combined_95_OR),
                        col = circlize::colorRamp2(c(-6, 0, 6), c("blue", "white", "#CC3333")),
                        cluster_rows = FALSE,
                        cluster_columns = FALSE,
                        cell_fun = function(j, i, x, y, width, height, fill) {
                          grid.text(combined_95_OR_char[i, j], x, y, gp = gpar(fontsize = 10))
                        },
                        column_names_rot = 0,
                        name = "log2(Odd's\nratio)")



combined_99 <- combined[which(combined$identity == '99'), ]

combined_99_OR <- reshape2::dcast(data = combined_99, formula = tax_level ~ category, value.var = "log2OR", fill = NA)
combined_99_BH <- reshape2::dcast(data = combined_99, formula = tax_level ~ category, value.var = "BH", fill = 1)

rownames(combined_99_OR) <- combined_99_OR$tax_level
combined_99_OR <- combined_99_OR[, -1]

rownames(combined_99_BH) <- combined_99_BH$tax_level
combined_99_BH <- combined_99_BH[, -1]

combined_99_OR_char <- apply(combined_99_OR, 2, function(x) { gsub(" *", "", sprintf("%.2f", x)) })
combined_99_OR_char <- apply(combined_99_OR_char, 2, function(x) { gsub(" *", "", x) })
combined_99_OR_char <- apply(combined_99_OR_char, 2, function(x) { gsub("NA", "", x) })
combined_99_OR_char[is.na(combined_99_OR_char)] <- ''

combined_99_OR[combined_99_BH >= 0.05] <- NA


ComplexHeatmap::Heatmap(matrix = as.matrix(combined_99_OR),
                        col = circlize::colorRamp2(c(-6, 0, 6), c("blue", "white", "#CC3333")),
                        cluster_rows = FALSE,
                        cluster_columns = FALSE,
                        cell_fun = function(j, i, x, y, width, height, fill) {
                          grid.text(combined_99_OR_char[i, j], x, y, gp = gpar(fontsize = 10))
                        },
                        column_names_rot = 0,
                        name = "log2(Odd's\nratio)")

