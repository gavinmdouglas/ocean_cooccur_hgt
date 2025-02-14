rm(list = ls(all.names = TRUE))

library(ggplot2)
library(reshape2)

# Proportion of comparisons that were putative hits for the BLAST and cluster-based approaches.
blast_norm_hits <- read.table('/mfs/gdouglas/projects/ocean_hgt_zenodo/putative_hgt/blast/cross_level_tallies_norm.tsv.gz',
                              stringsAsFactors = FALSE, sep = '\t', header = TRUE)
blast_norm_hits <- blast_norm_hits[, c("Inter.level", grep(".norm", colnames(blast_norm_hits), value = TRUE))]
blast_norm_hits_long <- reshape2::melt(data = blast_norm_hits, id.vars = "Inter.level")
blast_norm_hits_long <- blast_norm_hits_long[-which(blast_norm_hits_long$Inter.level %in% c("Species", "Strain")), ]
blast_norm_hits_long$variable <- as.character(blast_norm_hits_long$variable)
blast_norm_hits_long$variable[which(blast_norm_hits_long$variable == "blast95_99.norm")] <- "Identity >= 95% and < 99%"
blast_norm_hits_long$variable[which(blast_norm_hits_long$variable == "blast99.norm")] <- "Identity >= 99%"
blast_norm_hits_long$variable[which(blast_norm_hits_long$variable == "both.norm")] <- "All hits (Identity >= 95%)"
blast_norm_hits_long$Approach <- 'BLASTn'


cluster_norm_hits <- read.table('/mfs/gdouglas/projects/ocean_hgt_zenodo/putative_hgt/cluster/cross_level_tallies_norm.tsv.gz',
                              stringsAsFactors = FALSE, sep = '\t', header = TRUE)
cluster_norm_hits <- cluster_norm_hits[, c("Level", grep(".norm", colnames(cluster_norm_hits), value = TRUE))]
colnames(cluster_norm_hits)[1] <- "Inter.level"
cluster_norm_hits_long <- reshape2::melt(data = cluster_norm_hits, id.vars = "Inter.level")
cluster_norm_hits_long$variable <- as.character(cluster_norm_hits_long$variable)
cluster_norm_hits_long$variable[which(cluster_norm_hits_long$variable == "cluster_95_99.norm")] <- "Identity >= 95% and < 99%"
cluster_norm_hits_long$variable[which(cluster_norm_hits_long$variable == "cluster_99.norm")] <- "Identity >= 99%"
cluster_norm_hits_long$variable[which(cluster_norm_hits_long$variable == "both.norm")] <- "All hits (Identity >= 95%)"
cluster_norm_hits_long$Approach <- 'CD-HIT'

norm_hits_combined <- rbind(blast_norm_hits_long, cluster_norm_hits_long)

tax_levels <- c("Genus", "Family","Order", "Class", "Phylum", "Domain")
norm_hits_combined$Inter.level <- factor(norm_hits_combined$Inter.level, levels = tax_levels)

norm_hits_combined$variable <- factor(norm_hits_combined$variable, levels = c("All hits (Identity >= 95%)", "Identity >= 95% and < 99%", "Identity >= 99%"))

colnames(norm_hits_combined)[which(colnames(norm_hits_combined) == "variable")] <- "identity_cutoff"
colnames(norm_hits_combined)[which(colnames(norm_hits_combined) == "value")] <- "num_hgt_by_num_compare"


mixed_color <- (col2rgb("orange") +  col2rgb("cornflowerblue")) / 2
mixed_color_hex <- rgb(mixed_color[1,]/255, mixed_color[2,]/255, mixed_color[3,]/255)

norm_hits_plot <- ggplot(data = norm_hits_combined,
                         aes(x = Inter.level, y = num_hgt_by_num_compare, shape = Approach, colour = identity_cutoff)) +
  geom_point(size = 8, alpha=1) +
  scale_y_continuous(trans = "log10", labels = scales::label_number(), limits=c(0.0000005, 1)) +
  theme_bw() +
  scale_colour_manual(values = c(mixed_color_hex, "orange", "cornflowerblue")) +
  labs(colour = "Identity cut-off") +
  xlab("Inter-taxonomic level") +
  ylab(expression(atop(displaystyle(frac("No. putative HGT events", "No. genome comparisons")),
                       displaystyle("(log"[10] * " scale)")))) +
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5))

ggsave(plot = norm_hits_plot,
       filename = "/mfs/gdouglas/scripts/ocean_mag_hgt/display_items/blast_and_cluster_hits_norm.pdf",
       device = "pdf", width = 10, height = 6, units = "in", dpi=600)
