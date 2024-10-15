rm(list = ls(all.names = TRUE))

# Identify samples that have yet to be analyzed.
# And figure out the filter sizes of the ones analyzed so far.

library("readxl")

TableS1 <- read_excel("/mfs/gdouglas/projects/ocean_mags/metadata/OceanDNA_supp_metadata/Supp_File_S1.xlsx", skip=2)

Nico_subset <- read_excel("/mfs/gdouglas/projects/ocean_mags/metadata/OceanDNA_supp_metadata/Samples_to_keep.xlsx", skip=1)

TableS1_remaining <- TableS1[which(! TableS1$sample_name %in% Nico_subset$sample_name), ]

TableS1_remaining <- TableS1_remaining[which(TableS1_remaining$sample_type == "water"), ]

PRJEB6608_metadata <- read.table("/mfs/gdouglas/projects/ocean_mags/metadata/PRJEB6608_metadata.csv",
                                 header = TRUE, sep = ',', stringsAsFactors = FALSE)

PRJEB9740_metadata <- read.table("/mfs/gdouglas/projects/ocean_mags/metadata/PRJEB9740_metadata.csv",
                                 header = TRUE, sep = ',', stringsAsFactors = FALSE)

PRJEB1787_metadata <- read.table("/mfs/gdouglas/projects/ocean_mags/metadata/PRJEB1787_metadata.csv",
                                 header = TRUE, sep = ',', stringsAsFactors = FALSE)

PRJEB6608_unique_samples <- unique(PRJEB6608_metadata$Sample)
length(which(PRJEB6608_unique_samples %in% TableS1_remaining$sample_name))
length(PRJEB6608_unique_samples)

PRJEB9740_unique_samples <- unique(PRJEB9740_metadata$Sample)
length(which(PRJEB9740_unique_samples %in% TableS1_remaining$sample_name))
length(PRJEB9740_unique_samples)
PRJEB9740_unique_samples[which(! PRJEB9740_unique_samples %in% TableS1_remaining$sample_name)]

PRJEB1787_unique_samples <- unique(PRJEB1787_metadata$Sample)
length(which(PRJEB1787_unique_samples %in% TableS1_remaining$sample_name))
length(PRJEB1787_unique_samples)
PRJEB1787_unique_samples[which(! PRJEB1787_unique_samples %in% TableS1_remaining$sample_name)]

PRJEB6608_PRJEB1787_intersect_samples <- intersect(PRJEB6608_unique_samples, PRJEB1787_unique_samples)
length(which(PRJEB6608_PRJEB1787_intersect_samples %in% TableS1_remaining$sample_name))
length(PRJEB6608_PRJEB1787_intersect_samples)
PRJEB6608_PRJEB1787_intersect_samples[which(! PRJEB6608_PRJEB1787_intersect_samples %in% TableS1_remaining$sample_name)]

TableS1_remaining_nonTara <- TableS1_remaining
TableS1_remaining_nonTara <- TableS1_remaining_nonTara[-which(TableS1_remaining_nonTara$sample_name %in% PRJEB9740_unique_samples), ]
TableS1_remaining_nonTara <- TableS1_remaining_nonTara[-which(TableS1_remaining_nonTara$sample_name %in% PRJEB1787_unique_samples), ]
TableS1_remaining_nonTara <- as.data.frame(TableS1_remaining_nonTara)

TableS1_remaining_nonTara <- TableS1_remaining_nonTara[-grep("ERS492814", TableS1_remaining_nonTara$sample_name), ]

# Ignore remaining biofilm samples classified as "water".
TableS1_remaining_nonTara <- TableS1_remaining_nonTara[-which(TableS1_remaining_nonTara$division == "biofilm"), ]

write.table(x = TableS1_remaining_nonTara,
            file = "/mfs/gdouglas/projects/ocean_mags/metadata/OceanDNA_supp_metadata/OceanDNA_not_downloaded.tsv",
            sep = '\t', quote = FALSE, col.names = TRUE, row.names = FALSE)

write.table(x = TableS1_remaining_nonTara$sra_run,
            file = '/mfs/gdouglas/projects/ocean_mags/additional_round2/additional_run_ids/OceanDNA_nonTara_run_ids.txt',
            sep = '', quote = FALSE, col.names = FALSE, row.names = FALSE)

# Then to get the filter ranges for ALL water samples.

TableS1_water <- TableS1[which(TableS1$sample_type == "water"), ]
TableS1_water <- TableS1_water[-which(TableS1_water$division == "biofilm"), ]

filter <- data.frame(lower_filter=TableS1_water$lower_filter, upper_filter=TableS1_water$upper_filter)

write.table(x = filter,
            file = "/mfs/gdouglas/tmp/filter_cutoffs.txt",
            sep = '\t', quote = FALSE, col.names = TRUE, row.names = FALSE)


# Then also figure out which FASTQs Nico downloaded from the TARA datasets (and re-download any that are missing or empty)
nico_PRJEB6608_runs <- gsub('_1_clean.fastq.gz', '', list.files('/mfs/nicot/PRJEB6608-fastq_ftp-20231101-1543/cleaned', pattern = '_1_clean.fastq.gz'))
which(duplicated(nico_PRJEB6608_runs))

setdiff(nico_PRJEB6608_runs, PRJEB6608_metadata$Run)
setdiff(PRJEB6608_metadata$Run, nico_PRJEB6608_runs)



nico_PRJEB1787_runs <- gsub('_1_clean.fastq.gz', '', list.files('/mfs/nicot/PRJEB1787-fastq_ftp-20231101-1507/cleaned', pattern = '_1_clean.fastq.gz'))
which(duplicated(nico_PRJEB1787_runs))

setdiff(nico_PRJEB1787_runs, PRJEB1787_metadata$Run)
setdiff(PRJEB1787_metadata$Run, nico_PRJEB1787_runs)



nico_PRJEB9740_runs <- gsub('_1_clean.fastq.gz', '', list.files('/mfs/nicot/PRJEB9740-fastq_ftp-20231101-1513/cleaned', pattern = '_1_clean.fastq.gz'))
which(duplicated(nico_PRJEB9740_runs))

setdiff(nico_PRJEB9740_runs, PRJEB9740_metadata$Run)
setdiff(PRJEB9740_metadata$Run, nico_PRJEB9740_runs)

