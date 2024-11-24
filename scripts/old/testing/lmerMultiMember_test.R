rm(list = ls(all.names = TRUE))

library(dplyr)
library(lmerMultiMember)

combined_file="/mfs/gdouglas/projects/ocean_mags/networks/filtersplit/rangerdtl/metaG_hyperg_cooccur.freeliv.combined.tsv.gz"
cooccur_approach='hyperg'
hgt_tally_col='ranger_hgt_tallies'
median_diff_file="/mfs/gdouglas/projects/ocean_mags/networks/combined_tables/freeliv_present_metadata_median_pairwise_diff.tsv.gz"
outfolder="/mfs/gdouglas/tmp"
num_cores=8
suffix='.freeliv.rds'
keep_lower_levels = TRUE


if (! cooccur_approach %in% c('hyperg', 'simple', 'propr')) { stop('Co-occur approach must be hyperg, simple, or propr.') }
if (length(grep(cooccur_approach, combined_file)) == 0) { stop('Co-occur approach not present in combined_file?') }

combined_info <- read.table(combined_file, header=TRUE, sep="\t", stringsAsFactors = FALSE, row.names = 1)

if (! keep_lower_levels) {
  combined_info <- combined_info[which(! combined_info$diff_tax_level %in% c('Species', 'Strain')), ]
}
combined_info <- combined_info[which(rowSums(is.na(combined_info)) == 0), ]

if (cooccur_approach == 'hyperg') {
  combined_info$cooccur_ratio <- combined_info$cooccur_obs / combined_info$cooccur_exp
  combined_info$cooccur <- 0
  combined_info$cooccur[which(combined_info$cooccur_BH < 0.05 & combined_info$cooccur_ratio > 1)] <- 1
} else if (cooccur_approach == 'simple') {
  combined_info <- combined_info[which(! is.na(combined_info$cooccur_simple_cooccur)), ]
  combined_info$cooccur <- bestNormalize::orderNorm(x = combined_info$cooccur_simple_cooccur)$x.t
} else if (cooccur_approach == 'propr') {
  combined_info <- combined_info[which(! is.na(combined_info$cooccur_asso)), ]
  combined_info$cooccur <- bestNormalize::orderNorm(x = combined_info$cooccur_asso)$x.t
}

if (cooccur_approach == 'simple') {
  combined_info$hgt <- bestNormalize::orderNorm(x = combined_info[, hgt_tally_col])$x.t
  glmm_family = "gaussian"
} else {
  combined_info$hgt <- 0
  combined_info$hgt[which(combined_info[, hgt_tally_col] > 0)] <- 1
  glmm_family = "binomial"
}

combined_info$tip_dist_orderedNorm <- bestNormalize::orderNorm(x = combined_info$tip_dist)$x.t

Wp <- weights_from_columns(combined_info[, c("taxon_i", "taxon_j")])

m1 <- glmer(hgt ~ cooccur + (1 | genome_pair),
            family = "binomial",
            memberships = list(genome_pair = Wp),
            data = combined_info)

m1_summary <- summary(m1)
m1_summary$coefficients


# Random effects
m1_ranefs <- broom.mixed::tidy(m1, effects = "ran_vals", conf.int = TRUE) %>%
  .[.$group == "genome_pair", ] %>%
  .[order(.$estimate, decreasing = TRUE), ]

# plot top 10
m1_ranefs[1:10, ] %>%
  ggplot(aes(x = estimate, y = factor(level, level = level))) +
  geom_point(size = 3, color = "coral1") +
  geom_linerange(aes(xmin = conf.low, xmax = conf.high),
                 size = 1, color = "coral1") +
  theme_bw() + xlab("") + ylab("") +
  ggtitle("Increase in log(odds ratio) associated with player") +
  theme(axis.ticks = element_blank(),
        axis.text.x = element_text(size = 12, angle = 60, hjust = 1.0),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size = 14),
        rect = element_rect(fill = "transparent")) +
  scale_y_discrete(limits = rev) +
  geom_vline(xintercept = 0, color = "coral2", size = 1.0)

