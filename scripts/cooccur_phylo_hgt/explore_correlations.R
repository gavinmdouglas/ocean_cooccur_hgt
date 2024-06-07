rm(list = ls(all.names = TRUE))

library(ggplot2)
library(ggbeeswarm)

summary_stats <- function(x) {
  n <- length(x)
  mean_val <- mean(x, na.rm = TRUE)
  median_val <- median(x, na.rm = TRUE)
  sd_val <- sd(x, na.rm = TRUE)
  se_val <- sd_val / sqrt(n)
  wilcox_out <- wilcox.test(x)
  out <- c(n, mean_val, median_val, sd_val, se_val, wilcox_out$statistic, wilcox_out$p.value)
  names(out) <- c("n", "mean", "median", "sd", "se", "Wilcox_V", "Wilcox_P")
  return(out)
}

read_corr_info <- function(pairwise_file, partial_file, in_title) {
  pairwise_in <- read.table(pairwise_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  pairwise_in$comparison <- paste(pairwise_in$X, pairwise_in$Y, sep="v")
  
  category_map <- list()
  category_map[["cooccur_ratiovboth_gene_count"]] <- "Cooccurence vs.\nHGT gene count"
  category_map[["cooccur_ratiovtip_dist"]] <- "Cooccurence vs.\ntip dist."
  category_map[["both_gene_countvtip_dist"]] <- "HGT gene count\nvs. tip dist"
  
  category_map[["cooccur_simple_cooccurvboth_gene_count"]] <- "Cooccurence vs.\nHGT gene count"
  category_map[["cooccur_simple_cooccurvtip_dist"]] <- "Cooccurence vs.\ntip dist."

  category_map[["cooccur_assovboth_gene_count"]] <- "Cooccurence vs.\nHGT gene count"
  category_map[["cooccur_assovtip_dist"]] <- "Cooccurence vs.\ntip dist."
  
  pairwise_in$comparison_clean <- sapply(pairwise_in$comparison, function(x) { category_map[[x]] })
  
  pairwise_in$sig <- "Non-sig."
  pairwise_in$sig[which(pairwise_in$p.unc < 0.05)] <- "P < 0.05"
  
  pairwise_boxplots <- ggplot(data = pairwise_in, aes(x = comparison_clean, y = r)) +
                              geom_quasirandom(aes(colour=sig), alpha = 0.5) +
                              geom_boxplot(outlier.shape = NA, alpha=0.7) +
                              theme_bw() +
                              theme(panel.grid.major.x = element_blank(),
                                    panel.grid.minor.x = element_blank(),
                                    plot.title = element_text(hjust = 0.5)) +
                              scale_colour_manual(values=c("grey80", "salmon")) +
                              ylab("Spearman's correlation coefficient") +
                              xlab("Pairwise comparison") +
                              labs(colour="Test result") +
                              ggtitle(in_title)
  
  pairwise_r_summary <- aggregate(data = pairwise_in, x = r ~ comparison_clean, FUN = summary_stats)
  pairwise_r_summary$Category <- in_title
  
  partial_in <- read.table(partial_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  partial_in$Category <- in_title
  
  return(list(boxplots=pairwise_boxplots,
              partial_tab=partial_in,
              r_summary=pairwise_r_summary))
}

hyperg_w_lower <- read_corr_info(pairwise_file = "/mfs/gdouglas/projects/ocean_mags/coverm/network_working/corr_out/metaG_hyperG.w_lower.pairwise.tsv",
                                 partial_file = "/mfs/gdouglas/projects/ocean_mags/coverm/network_working/corr_out/metaG_hyperG.w_lower.partial.tsv",
                                 in_title = "Hypergeometric (with strains and species)")

hyperg_wo_lower <- read_corr_info(pairwise_file = "/mfs/gdouglas/projects/ocean_mags/coverm/network_working/corr_out/metaG_hyperG.pairwise.tsv",
                                 partial_file = "/mfs/gdouglas/projects/ocean_mags/coverm/network_working/corr_out/metaG_hyperG.partial.tsv",
                                 in_title = "Hypergeometric (without strains and species)")

simple_w_lower <- read_corr_info(pairwise_file = "/mfs/gdouglas/projects/ocean_mags/coverm/network_working/corr_out/metaG_simple.w_lower.pairwise.tsv",
                                 partial_file = "/mfs/gdouglas/projects/ocean_mags/coverm/network_working/corr_out/metaG_simple.w_lower.partial.tsv",
                                 in_title = "Simple co-occur (with strains and species)")

simple_wo_lower <- read_corr_info(pairwise_file = "/mfs/gdouglas/projects/ocean_mags/coverm/network_working/corr_out/metaG_simple.pairwise.tsv",
                                  partial_file = "/mfs/gdouglas/projects/ocean_mags/coverm/network_working/corr_out/metaG_simple.partial.tsv",
                                  in_title = "Simple co-occur (without strains and species)")

propr_w_lower <- read_corr_info(pairwise_file = "/mfs/gdouglas/projects/ocean_mags/coverm/network_working/corr_out/metaG_propr.w_lower.pairwise.tsv",
                                 partial_file = "/mfs/gdouglas/projects/ocean_mags/coverm/network_working/corr_out/metaG_propr.w_lower.partial.tsv",
                                 in_title = "propr co-abun. (with strains and species)")

propr_wo_lower <- read_corr_info(pairwise_file = "/mfs/gdouglas/projects/ocean_mags/coverm/network_working/corr_out/metaG_propr.pairwise.tsv",
                                  partial_file = "/mfs/gdouglas/projects/ocean_mags/coverm/network_working/corr_out/metaG_propr.partial.tsv",
                                  in_title = "propr co-abun. (without strains and species)")

# All pairwise boxplots.
hyperg_w_lower$boxplots
hyperg_wo_lower$boxplots
simple_w_lower$boxplots
simple_wo_lower$boxplots
propr_w_lower$boxplots
propr_wo_lower$boxplots

# Combine partial dataframes.
partial_w_lower <- hyperg_w_lower$partial_tab
partial_w_lower <- rbind(partial_w_lower, simple_w_lower$partial_tab)
partial_w_lower <- rbind(partial_w_lower, propr_w_lower$partial_tab)

partial_wo_lower <- hyperg_wo_lower$partial_tab
partial_wo_lower <- rbind(partial_wo_lower, simple_wo_lower$partial_tab)
partial_wo_lower <- rbind(partial_wo_lower, propr_wo_lower$partial_tab)

partial_w_lower$Category <- gsub(" \\(with strains and species\\)", "", partial_w_lower$Category)
partial_wo_lower$Category <- gsub(" \\(without strains and species\\)", "", partial_wo_lower$Category)

partial_w_lower$sig <- "Non-sig."
partial_w_lower$sig[which(partial_w_lower$p.val < 0.05)] <- "P < 0.05"

partial_w_lower_boxplots <- ggplot(data = partial_w_lower, aes(x = Category, y = r)) +
  geom_quasirandom(aes(colour=sig), alpha = 0.5) +
  geom_boxplot(outlier.shape = NA, alpha=0.7) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  scale_colour_manual(values=c("grey80", "salmon")) +
  ylab("Spearman's correlation coefficient (partial)") +
  xlab("Co-occurrence approach") +
  labs(colour="Test result") +
  ggtitle("Co-occur vs. HGT (with lower taxa)")

partial_wo_lower$sig <- "Non-sig."
partial_wo_lower$sig[which(partial_wo_lower$p.val < 0.05)] <- "P < 0.05"
partial_wo_lower_boxplots <- ggplot(data = partial_wo_lower, aes(x = Category, y = r)) +
  geom_quasirandom(aes(colour=sig), alpha = 0.5) +
  geom_boxplot(outlier.shape = NA, alpha=0.7) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  scale_colour_manual(values=c("grey80", "salmon")) +
  ylab("Spearman's correlation coefficient (partial)") +
  xlab("Co-occurrence approach") +
  labs(colour="Test result") +
  ggtitle("Co-occur vs. HGT (without lower taxa)")
