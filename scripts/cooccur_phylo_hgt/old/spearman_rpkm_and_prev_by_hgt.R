rm(list = ls(all.names = TRUE))

library(ComplexHeatmap)

abun_hgt <- read.table("/mfs/gdouglas/projects/ocean_mags/coverm/combined_tables/metaG_rpkm_allsamples_mean_by_hgt_partners.tsv.gz",
                       header = TRUE, sep = "\t", stringsAsFactors = FALSE)

hgt_types <- unique(abun_hgt$Category)

tax_levels <- unique(abun_hgt$Level)

vars_of_interest <- c("RPKM", "RPKM_present", "Samples_present")
spearman_out_template <- data.frame(matrix(NA, nrow = 3, ncol = 6))
rownames(spearman_out_template) <- vars_of_interest
colnames(spearman_out_template) <- c("variable", "hgt_type", "tax_level", "rho", "p", "N")
spearman_out_template$variable <- vars_of_interest

raw_combined <- list()

for (hgt_type in hgt_types) {
  
  for (tax_level in tax_levels) {
 
    tab_subset <- abun_hgt[which(abun_hgt$Category == hgt_type & abun_hgt$Level == tax_level), ]

    spearman_out <- spearman_out_template
    spearman_out$hgt_type <- hgt_type
    spearman_out$tax_level <- tax_level
    spearman_out$N <- nrow(tab_subset)
    
    for (var_of_interest in vars_of_interest) {
      cor_out <- cor.test(tab_subset[, var_of_interest],
                          tab_subset[, "Num_HGT_pairs"])
      spearman_out[var_of_interest, c("rho", "p")] <- c(cor_out$estimate, cor_out$p.value)
    }
    
    raw_combined[[paste(hgt_type, tax_level)]] <- spearman_out
    
  }
  
}

all_combined <- do.call(rbind, raw_combined)

all_combined$tax_level <- factor(all_combined$tax_level, levels = c("Genus", "Family", "Order", "Class", "Phylum", "Domain"))

rpkm_spearman <- all_combined[which(all_combined$variable == "RPKM"), ]

rpkm_spearman_rho <- reshape2::dcast(data = rpkm_spearman, formula = hgt_type ~ tax_level, value.var = "rho")
rpkm_spearmanp <- reshape2::dcast(data = rpkm_spearman, formula = hgt_type ~ tax_level, value.var = "p")

rownames(rpkm_spearman_rho) <- rpkm_spearman_rho$hgt_type
rpkm_spearman_rho <- rpkm_spearman_rho[, -1]

rownames(rpkm_spearmanp) <- rpkm_spearmanp$hgt_type
rpkm_spearmanp <- rpkm_spearmanp[, -1]

rpkm_spearman_rho_char <- apply(rpkm_spearman_rho, 2, function(x) { sprintf("%.2f", x) })
rpkm_spearman_rho[rpkm_spearmanp >= 0.05] <- NA


ComplexHeatmap::Heatmap(matrix = as.matrix(rpkm_spearman_rho),
                        col = circlize::colorRamp2(c(0, 1), c("white", "#CC3333")),
                        cluster_rows = FALSE,
                        cluster_columns = FALSE,
                        cell_fun = function(j, i, x, y, width, height, fill) {
                          grid.text(rpkm_spearman_rho_char[i, j], x, y, gp = gpar(fontsize = 10))
                        },
                        name = "Spearman rho")



prev_spearman <- all_combined[which(all_combined$variable == "Samples_present"), ]

prev_spearman_rho <- reshape2::dcast(data = prev_spearman, formula = hgt_type ~ tax_level, value.var = "rho")
prev_spearmanp <- reshape2::dcast(data = prev_spearman, formula = hgt_type ~ tax_level, value.var = "p")

rownames(prev_spearman_rho) <- prev_spearman_rho$hgt_type
prev_spearman_rho <- prev_spearman_rho[, -1]

rownames(prev_spearmanp) <- prev_spearmanp$hgt_type
prev_spearmanp <- prev_spearmanp[, -1]

prev_spearman_rho_char <- apply(prev_spearman_rho, 2, function(x) { sprintf("%.2f", x) })
prev_spearman_rho[prev_spearmanp >= 0.05] <- NA


ComplexHeatmap::Heatmap(matrix = as.matrix(prev_spearman_rho),
                        col = circlize::colorRamp2(c(0, 1), c("white", "#CC3333")),
                        cluster_rows = FALSE,
                        cluster_columns = FALSE,
                        cell_fun = function(j, i, x, y, width, height, fill) {
                          grid.text(prev_spearman_rho_char[i, j], x, y, gp = gpar(fontsize = 10))
                        },
                        name = "Spearman rho")