rm(list = ls(all.names = TRUE))

# Create quick plots highlighting difference in sample composition.

freeliving_map <- read.table('/mfs/gdouglas/projects/ocean_mags/metadata/OceanDNA_supp_metadata/OceanDNA_metadata_water_freeliving.tsv',
                             header=TRUE, sep='\t', stringsAsFactors = FALSE)

# Re-code one funky sample name to match presence table.
freeliving_map[which(freeliving_map$sample_name == 'ERS492821_ERS492814'), 'sample_name'] <- 'ERS492814'

lessfiltered_map <- read.table('/mfs/gdouglas/projects/ocean_mags/metadata/OceanDNA_supp_metadata/OceanDNA_metadata_water_notparticledepleted.tsv',
                               header=TRUE, sep='\t', stringsAsFactors = FALSE)

freeliv_data <- as.data.frame(table(freeliving_map$division))
names(freeliv_data) <- c("Division", "Count")

ggplot(freeliv_data, aes(x = reorder(Division, -Count), y = Count)) +
  geom_bar(stat = "identity", fill = "#2C7BB6") +
  geom_text(aes(label = Count), vjust = -0.5) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank()
  ) +
  labs(x = "Division", y = "Count", title = "Distribution by Division") +
  ggtitle('Free-living-associated samples')

lessfilt_data <- as.data.frame(table(lessfiltered_map$division))
names(lessfilt_data) <- c("Division", "Count")

ggplot(lessfilt_data, aes(x = reorder(Division, -Count), y = Count)) +
  geom_bar(stat = "identity", fill = "#2C7BB6") +
  geom_text(aes(label = Count), vjust = -0.5) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank()
  ) +
  labs(x = "Division", y = "Count", title = "Distribution by Division") +
  ggtitle('Less-filtered-associated samples')
freeliving_map <- read.table('/mfs/gdouglas/projects/ocean_mags/metadata/OceanDNA_supp_metadata/OceanDNA_metadata_water_freeliving.tsv',
                             header=TRUE, sep='\t', stringsAsFactors = FALSE)

# Re-code one funky sample name to match presence table.
freeliving_map[which(freeliving_map$sample_name == 'ERS492821_ERS492814'), 'sample_name'] <- 'ERS492814'

lessfiltered_map <- read.table('/mfs/gdouglas/projects/ocean_mags/metadata/OceanDNA_supp_metadata/OceanDNA_metadata_water_notparticledepleted.tsv',
                               header=TRUE, sep='\t', stringsAsFactors = FALSE)

# Create quick plots highlighting difference in sample composition.
xfreeliv_data <- as.data.frame(table(freeliving_map$division))
names(freeliv_data) <- c("Division", "Count")

ggplot(freeliv_data, aes(x = reorder(Division, -Count), y = Count)) +
  geom_bar(stat = "identity", fill = "#2C7BB6") +
  geom_text(aes(label = Count), vjust = -0.5) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank()
  ) +
  labs(x = "Division", y = "Count", title = "Distribution by Division") +
  ggtitle('Free-living-associated samples')

lessfilt_data <- as.data.frame(table(lessfiltered_map$division))
names(lessfilt_data) <- c("Division", "Count")

ggplot(lessfilt_data, aes(x = reorder(Division, -Count), y = Count)) +
  geom_bar(stat = "identity", fill = "#2C7BB6") +
  geom_text(aes(label = Count), vjust = -0.5) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank()
  ) +
  labs(x = "Division", y = "Count", title = "Distribution by Division") +
  ggtitle('Less-filtered-associated samples')
