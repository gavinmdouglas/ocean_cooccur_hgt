rm(list = ls(all.names = TRUE))

library(ggplot2)
library(cowplot)

combined_alldata <- read.table('/mfs/gdouglas/projects/ocean_mags/networks/allsamples/clusterbased/cooccur_and_tipdist.combined_subset_for_plotting.tsv.gz',
                             header=TRUE, sep='\t', stringsAsFactors = FALSE)

combined_lower <- combined_alldata[which(combined_alldata$diff_tax_level %in% c('Species', 'Strain')), ]

combined_higher <- combined_alldata[which(! combined_alldata$diff_tax_level %in% c('Species', 'Strain')), ]

hyperg_alldata <- ggplot(data = combined_alldata, aes(x = log2(hyperg_ratio), y = tip_dist)) +
  geom_hex(bins=20) +
  scale_fill_viridis_c() +
  theme_bw() +
  ggtitle('HyperG - all data') +
  ylab('Inter-tip distance')

hyperg_lower <- ggplot(data = combined_lower, aes(x = log2(hyperg_ratio), y = tip_dist)) +
  geom_hex(bins=20) +
  scale_fill_viridis_c() +
  theme_bw() +
  ggtitle('HyperG - species and strains') +
  ylab('Inter-tip distance')

hyperg_higher <- ggplot(data = combined_higher, aes(x = log2(hyperg_ratio), y = tip_dist)) +
  geom_hex(bins=20) +
  scale_fill_viridis_c() +
  theme_bw() +
  ggtitle('HyperG - genus and above') +
  ylab('Inter-tip distance')


simple_alldata <- ggplot(data = combined_alldata, aes(x = simple_cooccur, y = tip_dist)) +
  geom_hex(bins=20) +
  scale_fill_viridis_c() +
  theme_bw() +
  ggtitle('Simple co-occur - all data') +
  ylab('Inter-tip distance')

simple_lower <- ggplot(data = combined_lower, aes(x = simple_cooccur, y = tip_dist)) +
  geom_hex(bins=20) +
  scale_fill_viridis_c() +
  theme_bw() +
  ggtitle('Simple co-occur - species and strains') +
  ylab('Inter-tip distance')

simple_higher <- ggplot(data = combined_higher, aes(x = simple_cooccur, y = tip_dist)) +
  geom_hex(bins=20) +
  scale_fill_viridis_c() +
  theme_bw() +
  ggtitle('Simple co-occur - genus and above') +
  ylab('Inter-tip distance')





propr_alldata <- ggplot(data = combined_alldata, aes(x = propr_adja, y = tip_dist)) +
  geom_hex(bins=20) +
  scale_fill_viridis_c() +
  theme_bw() +
  ggtitle('propr - all data') +
  ylab('Inter-tip distance')

propr_lower <- ggplot(data = combined_lower, aes(x =propr_adja, y = tip_dist)) +
  geom_hex(bins=20) +
  scale_fill_viridis_c() +
  theme_bw() +
  ggtitle('propr - species and strains') +
  ylab('Inter-tip distance')

propr_higher <- ggplot(data = combined_higher, aes(x = propr_adja, y = tip_dist)) +
  geom_hex(bins=20) +
  scale_fill_viridis_c() +
  theme_bw() +
  ggtitle('propr - genus and above') +
  ylab('Inter-tip distance')


top_row <- plot_grid(hyperg_alldata, simple_alldata, propr_alldata, nrow = 1)
middle_row <- plot_grid(hyperg_lower, simple_lower, propr_lower, nrow = 1)
bottom_row <- plot_grid(hyperg_higher, simple_higher, propr_higher, nrow = 1)

plot_grid(top_row, middle_row, bottom_row, nrow = 3)
