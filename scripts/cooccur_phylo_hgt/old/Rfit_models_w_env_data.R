rm(list = ls(all.names = TRUE))

library(Rfit)
vars_of_interest <- c("depth", "latitude", "longitude", "temperature", "oxygen",
                      "tip_dist", "cooccur_simple_cooccur")

simple_cooccur_folder <- "/mfs/gdouglas/projects/ocean_mags/coverm/network_working_allsamples/metaG_simple_cooccur_prepped_tabs_with_env_mean_diffs/"
simple_cooccur_files <- list.files(simple_cooccur_folder, full.names = TRUE)

RAW_simple_cooccur_coef <- list()
RAW_simple_cooccur_R2 <- list()

for (infile in simple_cooccur_files) {

  taxon <- gsub(".tsv$", "", basename(infile))
  taxon_tab <- read.table(infile, header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = 1)

  hgt_var <- var(taxon_tab$both_gene_count)
  cooccur_var <- var(taxon_tab$cooccur_simple_cooccur)
  if (hgt_var == 0 | cooccur_var == 0) { next }

  tryCatch({

    model_out_wo_cooccur <- summary(rfit(formula = both_gene_count ~ depth + latitude + longitude + temperature + oxygen + tip_dist,
                                         data = taxon_tab))
    
    model_out_w_cooccur <- summary(rfit(formula = both_gene_count ~ depth + latitude + longitude + temperature + oxygen + tip_dist + cooccur_simple_cooccur,
                                        data = taxon_tab))
    
    RAW_simple_cooccur_R2[[taxon]] <- data.frame(taxon=taxon,
                                                 n=nrow(taxon_tab),
                                                 hgt_var=hgt_var,
                                                 cooccur_var=cooccur_var,
                                                 num_unique_cooccur=length(unique(taxon_tab$cooccur_simple_cooccur)),
                                                 nococcur=model_out_wo_cooccur$R2,
                                                 cooccur=model_out_w_cooccur$R2)
    
    model_out_wo_cooccur_coef <- data.frame(model_out_wo_cooccur$coefficients)
    model_out_w_cooccur_coef <- data.frame(model_out_w_cooccur$coefficients)
    model_out_wo_cooccur_coef$modeltype <- "nocooccur"
    model_out_w_cooccur_coef$modeltype <- "cooccur"
    
    RAW_simple_cooccur_coef[[taxon]] <- rbind(model_out_wo_cooccur_coef,
                                              model_out_w_cooccur_coef)
    
    RAW_simple_cooccur_coef[[taxon]]$taxon <- taxon
    
  }, error = function(e) {
    message("Error, skipping: ", taxon)
    next
  })
  
}

