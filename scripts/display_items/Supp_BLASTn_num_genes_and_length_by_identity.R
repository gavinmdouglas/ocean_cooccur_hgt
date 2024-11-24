rm(list = ls(all.names = TRUE))

library(ggplot2)
library(cowplot)
library(ggbeeswarm)

hit_info <- read.table(file = "/mfs/gdouglas/projects/ocean_hgt_zenodo/putative_hgt/blast/hit_gene_counts_and_lengths.tsv.gz",
                       header = TRUE, sep = "\t", stringsAsFactors = FALSE)

hit_info$Level <- factor(hit_info$Level, levels = c("Domain", "Phylum", "Class", "Order", "Family", "Genus"))

hit_info$Identity <- gsub("Identity ", "", hit_info$Identity)
hit_info$Identity <- factor(hit_info$Identity, levels=rev(unique(hit_info$Identity)))

gene_wilcox <- wilcox.test(hit_info$Num_covered_genes[which(hit_info$Identity == ">= 95% and < 99%")],
                           hit_info$Num_covered_genes[which(hit_info$Identity == ">= 99%")])

round(mean(hit_info$Num_covered_genes[which(hit_info$Identity == ">= 95% and < 99%")]), 1)
round(sd(hit_info$Num_covered_genes[which(hit_info$Identity == ">= 95% and < 99%")]), 1)

round(mean(hit_info$Num_covered_genes[which(hit_info$Identity == ">= 99%")]), 1)
round(sd(hit_info$Num_covered_genes[which(hit_info$Identity == ">= 99%")]), 1)

gene_wilcox_text <- paste0("Overall Wilcoxon Test\nW = ", format(gene_wilcox$statistic, trim=TRUE, big.mark=","),
                           "      \n", "P < 0.00001              ")

num_genes_plot <- ggplot(data = hit_info, aes(x = Num_covered_genes, y = Level, fill = Identity)) +
                    geom_quasirandom(aes(colour=Identity), alpha = 0.5, dodge.width=0.77, orientation="y") +
                    geom_boxplot(outlier.shape = NA, alpha = 0.5) +
                    theme_bw() +
                    xlab("Number of genes at least 90% covered by hit region") +
                    ylab("Inter-taxon level") +
                    scale_fill_manual(values=c("orange", "cornflowerblue")) +
                    scale_colour_manual(values=c("grey70", "grey70"), guide="none") +
                    guides(fill = guide_legend(reverse = TRUE)) +
                    annotate("text", x = Inf, y = -Inf,
                             label = gene_wilcox_text,
                             hjust = 1.1, vjust = -0.6,
                             size = 3)

length_wilcox <- wilcox.test(hit_info$length[which(hit_info$Identity == ">= 95% and < 99%")],
                             hit_info$length[which(hit_info$Identity == ">= 99%")])

length_wilcox_text <- paste0("Overall Wilcoxon Test\nW = ", format(length_wilcox$statistic, trim=TRUE, big.mark=","),
                           "      \n", "P < 0.00001              ")

length_plot <- ggplot(data = hit_info, aes(x = length, y = Level, fill = Identity)) +
                geom_quasirandom(aes(colour=Identity), alpha = 0.5, dodge.width=0.77, orientation="y") +
                geom_boxplot(outlier.shape = NA, alpha = 0.5) +
                theme_bw() +
                xlab("Hit region length") +
                ylab("Inter-taxon level") +
                scale_fill_manual(values=c("orange", "cornflowerblue")) +
                scale_colour_manual(values=c("grey70", "grey70"), guide="none") +
                guides(fill = guide_legend(reverse = TRUE)) +
                annotate("text", x = Inf, y = -Inf,
                         label = length_wilcox_text,
                         hjust = 1.1, vjust = -0.6,
                         size = 3)


combined_plot <- plot_grid(num_genes_plot, length_plot, labels=c('a', 'b'), nrow=2)

ggsave(plot = combined_plot,
       filename = "/mfs/gdouglas/scripts/ocean_mag_hgt/display_items/Supp_BLASTn_num_genes_and_length_by_identity.pdf",
       device = "pdf", width = 6, height = 7, units = "in", dpi=600)
