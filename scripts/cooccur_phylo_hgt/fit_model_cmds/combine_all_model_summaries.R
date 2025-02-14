rm(list = ls(all.names = TRUE))

all_rds <- list.files(path = '/mfs/gdouglas/projects/ocean_mags/glmm_working/',
                      pattern = '.rds$', recursive = TRUE, full.names = TRUE)

model_info <- list()
for (rdsfile in all_rds) {
  rds_in <- readRDS(rdsfile)
  model_info[[rdsfile]] <- list()
  model_info[[rdsfile]]$coef <- rds_in$coefficients
  if ('aic' %in% names(rds_in)) {
    model_info[[rdsfile]]$aic <- rds_in$aic
  } else {
    model_info[[rdsfile]]$aic <- as.numeric(rds_in$AICtab[1])
  }
  file_split <- strsplit(rdsfile, '/')[[1]]
  model_info[[rdsfile]]$genomes <- file_split[8]
  model_info[[rdsfile]]$tool_and_subset <- file_split[9]
  model_info[[rdsfile]]$approach <- file_split[10]
}

saveRDS(object = model_info, file = '/mfs/gdouglas/projects/ocean_mags/glmm_working/combined_summary.rds')


all_vif <- list.files(path = '/mfs/gdouglas/projects/ocean_mags/glmm_working/',
                      pattern = '.VIF.tsv$', recursive = TRUE, full.names = TRUE)

vif_info <- list()
for (filepath in all_vif) {
  vif_info[[filepath]] <- list()
  file_split <- strsplit(filepath, '/')[[1]]
  vif_info[[filepath]]$genomes <- file_split[8]
  vif_info[[filepath]]$tool_and_subset <- file_split[9]
  vif_info[[filepath]]$approach <- file_split[10]
  vif_info[[filepath]]$vif <- read.table(filepath, header=TRUE, sep='\t', stringsAsFactors = FALSE)
}

saveRDS(object = vif_info, file = '/mfs/gdouglas/projects/ocean_mags/glmm_working/all_vif.rds')

