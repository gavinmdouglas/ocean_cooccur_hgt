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
  category_map[["cooccur_ratiovranger_hgt_tallies"]] <- "Cooccurence vs.\nHGT gene count"
  category_map[["cooccur_ratiovtip_dist"]] <- "Cooccurence vs.\ntip dist."
  category_map[["ranger_hgt_talliesvtip_dist"]] <- "HGT gene count\nvs. tip dist"
  
  category_map[["cooccur_simple_cooccurvranger_hgt_tallies"]] <- "Cooccurence vs.\nHGT gene count"
  category_map[["cooccur_simple_cooccurvtip_dist"]] <- "Cooccurence vs.\ntip dist."

  category_map[["cooccur_assovranger_hgt_tallies"]] <- "Cooccurence vs.\nHGT gene count"
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


hyperg_lower <- read_corr_info(pairwise_file = "/mfs/gdouglas/projects/ocean_mags/coverm/network_working/ranger_hgt/corr_out/metaG_hyperG.pairwise.tsv",
                                 partial_file = "/mfs/gdouglas/projects/ocean_mags/coverm/network_working/ranger_hgt/corr_out/metaG_hyperG.partial.tsv",
                                 in_title = "Hypergeometric (RANGER-DTL results)")

simple_lower <- read_corr_info(pairwise_file = "/mfs/gdouglas/projects/ocean_mags/coverm/network_working/ranger_hgt/corr_out/metaG_simple.pairwise.tsv",
                                  partial_file = "/mfs/gdouglas/projects/ocean_mags/coverm/network_working/ranger_hgt/corr_out/metaG_simple.partial.tsv",
                                  in_title = "Simple co-occur (RANGER-DTL results)")


propr_lower <- read_corr_info(pairwise_file = "/mfs/gdouglas/projects/ocean_mags/coverm/network_working/ranger_hgt/corr_out/metaG_propr.pairwise.tsv",
                                  partial_file = "/mfs/gdouglas/projects/ocean_mags/coverm/network_working/ranger_hgt/corr_out/metaG_propr.partial.tsv",
                                  in_title = "propr co-abun. (RANGER-DTL results)")

# All pairwise boxplots.
hyperg_lower$boxplots
simple_lower$boxplots
propr_lower$boxplots


# Combine partial dataframes.
partial_combined <- rbind(hyperg_lower$partial_tab, simple_lower$partial_tab)
partial_combined <- rbind(partial_combined, propr_lower$partial_tab)

partial_combined$sig <- "Non-sig."
partial_combined$sig[which(partial_combined$p.val < 0.05)] <- "P < 0.05"

partial_combined$Category <- gsub(" \\(RANGER-DTL results\\)", "", partial_combined$Category)

partial_combined_boxplots <- ggplot(data = partial_combined, aes(x = Category, y = r)) +
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
  ggtitle("Co-occur vs. HGT (RANGER-DTL results)")


partial_combined_by_species <- aggregate(data = partial_combined, x =  r ~ species + Category, FUN=mean)

partial_combined_by_species_boxplots <- ggplot(data = partial_combined_by_species, aes(x = Category, y = r)) +
  geom_quasirandom(alpha = 0.5) +
  geom_boxplot(outlier.shape = NA, alpha=0.7) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  scale_colour_manual(values=c("grey80", "salmon")) +
  ylab("Spearman's correlation coefficient (partial)") +
  xlab("Co-occurrence approach") +
  labs(colour="Test result") +
  ggtitle("Co-occur vs. HGT (RANGER-DTL results, by species)")
