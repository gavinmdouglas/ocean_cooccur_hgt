rm(list = ls(all.names = TRUE))

library(ggplot2)
library(reshape2)

# Basic table and visualization of % comparisons that were hits by taxonomic level and identity cut-off.
num_comparisons <- read.table('/mfs/gdouglas/projects/ocean_mags/ms_summary_data/overall_summary/num_comparisons_per_inter.level.tsv.gz',
                              header = TRUE, sep = '\t', stringsAsFactors = FALSE, row.names = 1)

hgt_tab <- read.table('/mfs/gdouglas/projects/ocean_mags/ms_summary_data/clusterbased_summary/all_best_hits.tsv.gz',
                      header = TRUE, sep = '\t', stringsAsFactors = FALSE)

hgt_tab_95 <- hgt_tab[which(hgt_tab$identity >= 95.0 & hgt_tab$identity < 99.0), ]
hgt_tab_99 <- hgt_tab[which(hgt_tab$identity >= 99.0), ]



num_hits <- data.frame(Inter.level = tally_tab$Level[1:8],
                       num_comparisons = num_comparisons$Num_comparisons,
                       blast95_99=tally_tab[which(tally_tab$Identity == "Identity >= 95% and < 99%"), "Num_hits"],
                       blast99=tally_tab[which(tally_tab$Identity == "Identity >= 99%"), "Num_hits"])

num_hits$both <- num_hits$blast95_99 + num_hits$blast99

num_hits$blast95_99.norm <- num_hits$blast95_99 / num_comparisons$Num_comparisons
num_hits$blast99.norm <- num_hits$blast99 / num_comparisons$Num_comparisons
num_hits$both.norm <- num_hits$both / num_comparisons$Num_comparisons

write.table(x = num_hits,
            file = "~/scripts/ocean_mag_hgt/data/blast_summary/cross_level_tallies_norm.tsv",
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

norm_hits <- num_hits[, c("Inter.level", grep(".norm", colnames(num_hits), value = TRUE))]

norm_hits_long <- reshape2::melt(data = norm_hits, id.vars = "Inter.level")

norm_hits_long <- norm_hits_long[-which(norm_hits_long$Inter.level %in% c("Species", "Strain")), ]

norm_hits_long$Inter.level <- factor(norm_hits_long$Inter.level, levels = rev(num_hits$Inter.level))

norm_hits_long$variable <- as.character(norm_hits_long$variable)
norm_hits_long$variable[which(norm_hits_long$variable == "blast95_99.norm")] <- "Identity >= 95% and < 99%"
norm_hits_long$variable[which(norm_hits_long$variable == "blast99.norm")] <- "Identity >= 99%"
norm_hits_long$variable[which(norm_hits_long$variable == "both.norm")] <- "All hits (Identity >= 95%)"
norm_hits_long$variable <- factor(norm_hits_long$variable, levels = c("All hits (Identity >= 95%)", "Identity >= 95% and < 99%", "Identity >= 99%"))

norm_hits_plot <- ggplot(data = norm_hits_long,
                         aes(x = Inter.level, y = value, colour = variable)) +
                    geom_point(size = 5) +
                    scale_y_continuous(trans = "log10", labels = scales::label_number(), limits=c(0.00001, 0.1)) +
                    theme_bw() +
                    scale_colour_manual(values = c("black", "orange", "cornflowerblue")) +
                    labs(colour = "BLAST cut-off") +
                    xlab("Inter-taxonomic level") +
                    ylab("HGT\nrate") +
                    theme(axis.title.y = element_text(angle = 0, vjust = 0.5))

ggsave(plot = norm_hits_plot,
       filename = "/mfs/gdouglas/scripts/ocean_mag_hgt/display_items/blast_hits_norm.pdf",
       device = "pdf", width = 8, height = 5, units = "in", dpi=600)
