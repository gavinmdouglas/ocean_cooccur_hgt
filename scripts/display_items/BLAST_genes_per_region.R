
### OLD CODE THAT CAN BE USED TO GET THIS RESULT!

# Number of gene hits per pairwise BLAST region {.tabset}

Clearly some strains are just highly similar over huge regions, unsurprisingly.

Perhaps the biggest point of interest is that the average number of genes per region is actually pretty stable at the genus level and higher, but that regions that are more divergent tend to be smaller (which would make sense if they were older).

## With outliers

```{r genes_per_region_w_outliers}
hit_counts_per_region <- readRDS("/mfs/gdouglas/projects/ocean_mags/summary_files/per_hit_gene_counts.rds")

hit_counts_per_region <- hit_counts_per_region[which(! hit_counts_per_region$Level %in% c('Species', 'Strain')), ]
hit_counts_per_region$Level <- factor(hit_counts_per_region$Level, levels = c("Domain", "Phylum", "Class", "Order", "Family", "Genus"))

#hit_counts_per_region$Level <- factor(hit_counts_per_region$Level, levels = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain"))

ggplot(data = hit_counts_per_region, aes(x = Num_covered_genes, y = Level)) +
  geom_boxplot() +
  theme_bw() +
  facet_wrap(. ~ Identity) +
  xlab("Number of genes at least 90% covered by hit region\n(Average of genes in overlapping region of both genomes)") +
  ylab("Inter-taxon level")
```


## Outliers removed

```{r genes_per_region_no_outliers}
ggplot(data = hit_counts_per_region, aes(x = Num_covered_genes, y = Level)) +
  geom_boxplot(outlier.shape = NA) +
  theme_bw() +
  facet_wrap(. ~ Identity) +
  xlab("Number of genes at least 90% covered by hit region\n(Average of genes in overlapping region of both genomes)") +
  ylab("Inter-taxon level") +
  scale_x_continuous(limits = c(0, 20))
```

## Outliers removed -- reformatted
```{r genes_per_region_no_outliers_reformat}
ggplot(data = hit_counts_per_region, aes(x = Num_covered_genes, y = Level, fill = Identity)) +
  geom_boxplot(outlier.shape = NA) +
  theme_bw() +
  xlab("Number of genes at least 90% covered by hit region\n(Average of genes in overlapping region of both genomes)") +
  ylab("Inter-taxon level") +
  scale_fill_manual(values=c("green", "cornflowerblue")) +
  scale_x_continuous(limits = c(0, 20))
```

## Wilcoxon tests
```{r genes_wilcox_compare}

wilcox_summary <- data.frame(matrix(NA, nrow = num_tax_levels, ncol = 4))
rownames(wilcox_summary) <- c("Genus", "Family", "Order", "Class", "Phylum", "Domain")
#rownames(wilcox_summary) <- c("Strain", "Species", "Genus", "Family", "Order", "Class", "Phylum", "Domain")
colnames(wilcox_summary) <- c("mean95", "mean99", "W", "P")

for (tax_level in unique(hit_counts_per_region$Level)) {
  subset99 <- hit_counts_per_region[which(hit_counts_per_region$Identity == "Identity >= 99%" & hit_counts_per_region$Level == tax_level), "Num_covered_genes"]
  subset95 <- hit_counts_per_region[which(hit_counts_per_region$Identity == "Identity >= 95% and < 99%" & hit_counts_per_region$Level == tax_level), "Num_covered_genes"]

  wilcox_out <- wilcox.test(subset95, subset99)

  wilcox_summary[tax_level, "mean95"] <- mean(subset95)
  wilcox_summary[tax_level, "mean99"] <- mean(subset99)
  wilcox_summary[tax_level, "W"] <- wilcox_out$statistic
  wilcox_summary[tax_level, "P"] <- wilcox_out$p.value
}

kable(wilcox_summary) %>%
  kable_styling(full_width = FALSE)
```