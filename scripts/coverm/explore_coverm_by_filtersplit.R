rm(list = ls(all.names = TRUE))

# See how sample metadata and CoverM presence/absence,
# as well as abundance, differs between the free-living and less-filtered metagenomics samples.

write_gzip_table_w_rowcol(x = coverm_lessfiltered_presence,
                          outfile = )

write_gzip_table_w_rowcol(x = coverm_lessfiltered_rpkm,
                          outfile = '~/projects/ocean_mags/networks/combined_tables/metaG_rpkm_lessfiltered.tsv.gz')

freeliv_presence <- read.table('~/projects/ocean_mags/networks/combined_tables/metaG_presence_freeliv.tsv.gz',
                               header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names=1)

freeliv_rpkm <- read.table('~/projects/ocean_mags/networks/combined_tables/metaG_rpkm_freeliv.tsv.gz',
                               header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names=1)

lessfiltered_presence <- read.table('~/projects/ocean_mags/networks/combined_tables/metaG_presence_lessfiltered.tsv.gz',
                               header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names=1)

lessfiltered_rpkm <- read.table('~/projects/ocean_mags/networks/combined_tables/metaG_rpkm_lessfiltered.tsv.gz',
                           header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names=1)


metadata <- read.table("~/projects/ocean_mags/metadata/OceanDNA_supp_metadata/Supp_File_S1_water_samples.tsv",
                       stringsAsFactors = FALSE, header = TRUE, sep = "\t")
rownames(metadata) <- metadata$sample_name

freeliv_metadata <- metadata[rownames(freeliv_presence), ]
lessfiltered_metadata <- metadata[rownames(lessfiltered_presence), ]

boxplot(as.numeric(freeliv_metadata$oxygen), as.numeric(lessfiltered_metadata$oxygen), names=c("Free-living", "Less-filtered"), ylab = "Oxygen")
boxplot(as.numeric(freeliv_metadata$temperature), as.numeric(lessfiltered_metadata$temperature), names=c("Free-living", "Less-filtered"), ylab = "Temperature")
boxplot(as.numeric(freeliv_metadata$salinity), as.numeric(lessfiltered_metadata$salinity), names=c("Free-living", "Less-filtered"), ylab = "Salinity")
boxplot(as.numeric(freeliv_metadata$longigute), as.numeric(lessfiltered_metadata$longigute), names=c("Free-living", "Less-filtered"), ylab = "Longitude")
boxplot(as.numeric(freeliv_metadata$latitude), as.numeric(lessfiltered_metadata$latitude), names=c("Free-living", "Less-filtered"), ylab = "Latitude")
boxplot(log10(as.numeric(freeliv_metadata$depth) + 1), log10(as.numeric(lessfiltered_metadata$depth) + 1), names=c("Free-living", "Less-filtered"), ylab = "log10(Depth + 1)")
