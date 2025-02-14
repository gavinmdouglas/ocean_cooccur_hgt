rm(list = ls(all.names = TRUE))

# Subset OceanDNA table to only water samples of interest.

library("readxl")

TableS1 <- read_excel("/mfs/gdouglas/projects/ocean_mags/metadata/OceanDNA_supp_metadata/Supp_File_S1.xlsx", skip=2)

TableS1 <- TableS1[which(TableS1$sample_type == "water"), ]
TableS1 <- TableS1[-which(TableS1$division == "biofilm"), ]

write.table(x = TableS1,
            file = '/mfs/gdouglas/projects/ocean_mags/metadata/OceanDNA_supp_metadata/Supp_File_S1_water_samples.tsv',
            sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE)
