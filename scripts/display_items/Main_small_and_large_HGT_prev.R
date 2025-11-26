rm(list = ls(all.names = TRUE))

library(ggplot2)
library(ggbeeswarm)
library(reshape2)

combined_summary <- read.table('/mfs/gdouglas/projects/ocean_hgt_zenodo/hgt_prev_analyses/Tara_fraction_sample_matched_HGT_prev.tsv.gz',
                                header=TRUE, sep = '\t', stringsAsFactors = FALSE)

combined_summary <- combined_summary[which(combined_summary$lower_total_genomes >= 10 & combined_summary$upper_total_genomes >= 10), ]

# Also decided to remove all "prop_mult" panels, it we decided this was too much information,
# and unclear how to interpret (at the moment at least!)

col_to_rm <- grep("filt_group|_sample|total_genomes|mean_genes|_prop_mult_", colnames(combined_summary), value = TRUE)
combined_summary <- combined_summary[, which(! colnames(combined_summary) %in% col_to_rm)]

lower_cols <- grep("^lower_", names(combined_summary), value = TRUE)
upper_cols <- grep("^upper_", names(combined_summary), value = TRUE)
var_names <- gsub("^lower_", "", lower_cols)

results_table <- data.frame(
  variable = character(),
  lower_mean = numeric(),
  upper_mean = numeric(),
  mean_diff = numeric(),
  wilcox_p = numeric(),
  stringsAsFactors = FALSE
)

for (i in 1:length(var_names)) {
  var_name <- var_names[i]
  lower_col <- paste0("lower_", var_name)
  upper_col <- paste0("upper_", var_name)

  lower_vals <- combined_summary[, lower_col]
  upper_vals <- combined_summary[, upper_col]

  lower_mean <- mean(lower_vals, na.rm = TRUE)
  upper_mean <- mean(upper_vals, na.rm = TRUE)
  mean_diff <- upper_mean - lower_mean

  # Run paired Wilcoxon test
  wilcox_result <- wilcox.test(lower_vals, upper_vals, paired = TRUE)
  wilcox_p <- wilcox_result$p.value

  results_table <- rbind(results_table, data.frame(
    variable = var_name,
    lower_mean = lower_mean,
    upper_mean = upper_mean,
    mean_diff = mean_diff,
    wilcox_p = wilcox_p
  ))

}

results_table$wilcox_bh <- p.adjust(results_table$wilcox_p, 'BH')

remap_variables <- list(prop_hgt = "Horizontal gene transfer call",
                        prop_proMGE = "proMGE hit",
                        prop_plasmid = "Gene plasmid-annotated",
                        prop_virus = "Gene virus-annotated",
                        prop_provirus = "Provirus-annotated scaffold")
var_names <- names(remap_variables)

sig_data_sig <- data.frame(variable = results_table$variable[results_table$wilcox_bh < 0.05],
                           significance = "*",
                           stringsAsFactors = FALSE)

sig_data_ns <- data.frame(variable = results_table$variable[results_table$wilcox_bh >= 0.05],
                          significance = "n.s.",
                          stringsAsFactors = FALSE)

sig_data_sig$variable_label <- unlist(remap_variables[sig_data_sig$variable])
sig_data_sig$variable_label <- factor(sig_data_sig$variable_label, levels = unlist(remap_variables))

sig_data_ns$variable_label <- unlist(remap_variables[sig_data_ns$variable])
sig_data_ns$variable_label <- factor(sig_data_ns$variable_label, levels = unlist(remap_variables))

plot_data_list <- list()

for (i in 1:length(var_names)) {
  var_name <- var_names[i]
  lower_col <- paste0("lower_", var_name)
  upper_col <- paste0("upper_", var_name)

  lower_vals <- combined_summary[, lower_col]
  upper_vals <- combined_summary[, upper_col]

  temp_data <- data.frame(
    value = c(lower_vals, upper_vals),
    group = c(rep("lower", length(lower_vals)), rep("upper", length(upper_vals))),
    variable = var_name,
    stringsAsFactors = FALSE
  )

  plot_data_list[[i]] <- temp_data
}

plot_data <- do.call(rbind, plot_data_list)

plot_data$variable_label <- unlist(remap_variables[plot_data$variable])
plot_data$variable_label <- factor(plot_data$variable_label, levels = unlist(remap_variables))

plot_data$group <- ifelse(plot_data$group == "lower", yes = "Free-\nliving", no = "Particle-\nattached")

plot_data$group <- factor(plot_data$group, levels=c('Free-\nliving', 'Particle-\nattached'))

all_boxplots <- ggplot(plot_data, aes(x = group, y = value)) +
                        geom_quasirandom(alpha = 0.6, width = 0.3) +
                        geom_boxplot(alpha = 0.3, outlier.shape = NA, width = 0.5) +
                        facet_wrap(~ variable_label, scales = "free_y", axes="all_x") +
                        geom_text(data = sig_data_sig, aes(x = 1.5, y = Inf, label = significance),
                                  vjust = 1.2, size = 6, color = "red", inherit.aes = FALSE) +
                        geom_text(data = sig_data_ns, aes(x = 1.5, y = Inf, label = significance),
                                  vjust = 1.2, size = 3, color = "black", inherit.aes = FALSE) +
                        labs(x = "Sample size fraction",
                             y = "Proportion of genomes with at least one instance") +
                        theme_bw() +
                        theme(text = element_text(colour="black"),
                              axis.text = element_text(colour="black"),
                              strip.background = element_blank(),
                              strip.text = element_text(colour = "black"))

ggsave(filename = "/mfs/gdouglas/scripts/ocean_cooccur_hgt/display_items/Main_Figure7.pdf",
       plot = all_boxplots, height=4, width=6, units = "in", device="pdf", dpi=600)
