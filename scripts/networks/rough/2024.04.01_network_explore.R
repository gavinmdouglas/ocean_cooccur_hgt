rm(list = ls(all.names = TRUE))

library(reshape2)
library(parallel)
library(data.table)

bioprojects <- c('PRJEB1787', 'PRJEB97740', 'PRJEB6608')

bioproject_tallies <- data.frame(bioprojects = bioprojects,
                                 sample_number = NA)
rownames(bioproject_tallies) <- bioprojects

raw_in <- list()

for (bioproject in bioprojects) {
  
  all_files <- list.files(paste0("~/projects/water_mags/coverm/", bioproject), 
                          pattern = "cov.gz", full.names = TRUE)
  
  bioproject_tallies[bioproject, 'sample_number'] <- length(all_files)
  
  
  raw_in[[bioproject]] <- mclapply(all_files,
                     function(file) {
                       cov_in <- read.table(file, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
                       
                       colnames(cov_in) <- c('genome', 'mean', 'rpkm', 'trimmed_mean', 'relabun', 'covered_bases', 'length')
                       
                       cov_in$breadth <- cov_in$covered_bases / cov_in$length
                       
                       orig_col <- colnames(cov_in)
                       
                       mgs_id <- gsub('\\..*$', '', basename(file))
                       cov_in$mgs_sample <- mgs_id
                       
                       cov_in$bioproject <- bioproject
                       
                       return(cov_in[, c('bioproject', 'mgs_sample', orig_col)])

                      },
                     mc.cores=64)
}

# Split into metaG and metaT datasets.
combined_dna <- rbind(do.call(rbind, raw_in$PRJEB1787),
                      do.call(rbind, raw_in$PRJEB97740))
combined_rna <- do.call(rbind, raw_in$PRJEB6608)

# Then analyze data - starting with the DNA projects.
combined_dna$bioproject <- factor(combined_dna$bioproject, levels = bioprojects)

combined_dna$present <- 0

combined_dna[which(combined_dna$breadth >= 0.25), "present"] <- 1

dna_presence_matrix <- reshape2::dcast(data = combined_dna,
                                      formula = mgs_sample ~ genome,
                                      value.var = "present",
                                      fill = 0)

rownames(dna_presence_matrix) <- dna_presence_matrix$mgs_sample
dna_presence_matrix <- dna_presence_matrix[, which(colnames(dna_presence_matrix) != "mgs_sample")]

dna_sample_tallies <- rowSums(dna_presence_matrix)
dna_presence_matrix_filt <- dna_presence_matrix[which(dna_sample_tallies >= 50), ]

dna_taxon_presence_tallies <- colSums(dna_presence_matrix_filt)
dna_presence_matrix_filt <- dna_presence_matrix_filt[, which(dna_taxon_presence_tallies >= 5)]

dna_presence_matrix_filt_vecs <- lapply(dna_presence_matrix_filt, function(x) { which(x == 1) })

dna_presence_matrix_filt_vecs_lengths <- sapply(dna_presence_matrix_filt_vecs, length)

# Get index combinations to loop over.
dna_i_combinations <- combn(1:length(dna_presence_matrix_filt_vecs), 2)
dna_i_combinations_list <- mclapply(1:ncol(dna_i_combinations),
                                    function(i) { as.integer(dna_i_combinations[, i]) },
                                    mc.cores = 32)

combo_batch_ranges <- c(seq(from=1, to=length(dna_i_combinations_list), by=1000000), length(dna_i_combinations_list))

for (batch_i in 1:(length(combo_batch_ranges) - 1)) {
  
  start_index <- combo_batch_ranges[batch_i]
  end_index <- combo_batch_ranges[batch_i + 1] - 1
  
  if (batch_i == (length(combo_batch_ranges) - 1)) {
    end_index <- combo_batch_ranges[batch_i + 1]
  }
  
  dna_simple_cooccur_raw <- mclapply(dna_i_combinations_list[start_index:end_index],
                                     function(combo) {
                                       i = combo[1]
                                       j = combo[2]
                                       num_intersect = length(intersect(dna_presence_matrix_filt_vecs[[i]],
                                                                        dna_presence_matrix_filt_vecs[[j]]))
                                       
                                       min_num = min(dna_presence_matrix_filt_vecs_lengths[[i]],
                                                     dna_presence_matrix_filt_vecs_lengths[[j]])
                                       
                                       return(list(taxoni=i,
                                                   taxonj=j,
                                                   num_intersect=num_intersect,
                                                   min_num=min_num,
                                                   simple_cooccur=num_intersect/min_num))
                                     },
                                     mc.cores = 64)
  
  dna_simple_cooccur <- data.table::rbindlist(dna_simple_cooccur_raw)
  dna_simple_cooccur <- dna_simple_cooccur[which(dna_simple_cooccur$simple_cooccur > 0), ]
  dna_simple_cooccur$taxoni <- names(dna_presence_matrix_filt_vecs)[dna_simple_cooccur$taxoni]
  dna_simple_cooccur$taxonj <- names(dna_presence_matrix_filt_vecs)[dna_simple_cooccur$taxonj]
  
  saveRDS(object = dna_simple_cooccur,
          file = paste("/mfs/gdouglas/projects/water_mags/coverm/network_working/intermediate/dna_simple_cooccur_set", as.character(batch_i), ".rds", sep = ""))
  
}

output_rds_files <- list.files("/mfs/gdouglas/projects/water_mags/coverm/network_working/intermediate/", full.names = TRUE)

dna_simple_cooccur_combined <- readRDS(output_rds_files[1])

for (i in 2:length(output_rds_files)) {
  dna_simple_cooccur_combined <- rbind(dna_simple_cooccur_combined, readRDS(output_rds_files[i]))
}

saveRDS(object = dna_simple_cooccur_combined,
        file = "/mfs/gdouglas/projects/water_mags/coverm/network_working//dna_simple_cooccur.rds")







dna_simple_cooccur_raw <- mclapply(dna_i_combinations_list,
                                   function(combo) {
                                     i = combo[1]
                                     j = combo[2]
                                     num_intersect = length(intersect(dna_presence_matrix_filt_vecs[[i]],
                                                                      dna_presence_matrix_filt_vecs[[j]]))

                                     min_num = min(dna_presence_matrix_filt_vecs_lengths[[i]],
                                                   dna_presence_matrix_filt_vecs_lengths[[j]])
                                     
                                     return(list(taxoni=i,
                                                       taxonj=j,
                                                       num_intersect=num_intersect,
                                                       min_num=min_num,
                                                       simple_cooccur=num_intersect/min_num))
                                    },
                                   mc.cores = 64)

dna_simple_cooccur <- rbindlist(l = dna_simple_cooccur_raw, use.names = TRUE)

dna_simple_cooccur_df <- do.call(rbind,
                                 lapply(dna_simple_cooccur_raw, function(x) data.frame(t(x), stringsAsFactors = FALSE)))

dna_simple_cooccur <- rbindlist(l = dna_simple_cooccur_raw, use.names = TRUE)

plan(multisession, workers = 64)
dna_simple_cooccur <- future_rbindlist(dna_simple_cooccur_raw, use.names = TRUE, fill = TRUE)


dna_simple_cooccur <- do.call(rbind, dna_simple_cooccur_raw)

dna_simple_cooccur$taxon1 <- names(dna_presence_matrix_filt_vecs)[dna_simple_cooccur$]







dna_simple_cooccur_raw <- list()

for (combo_i in 1:length(dna_i_combinations_list)) {
  
  combo = dna_i_combinations_list[[combo_i]]
  i = combo[1]
  j = combo[2]
  num_intersect = length(intersect(dna_presence_matrix_filt_vecs[[i]],
                                   dna_presence_matrix_filt_vecs[[j]]))
  
  min_num = min(dna_presence_matrix_filt_vecs_lengths[[i]],
                dna_presence_matrix_filt_vecs_lengths[[j]])
  
  dna_simple_cooccur_raw[[combo_i]] <- data.frame(taxoni=i,
                                                  taxonj=j,
                                                  num_intersect=num_intersect,
                                                  min_num=min_num,
                                                  simple_cooccur=num_intersect/min_num)
  
}













# Try co-occur score.
library(cooccur)


test <- as.matrix(presence_matrix_filt)

tmp <- test[, 1:10]
cooccur_out <- cooccur(mat=t(tmp),
                       type="spp_site",
                       thresh=TRUE,
                       spp_names=TRUE)

cooccur_out <- cooccur_parallel(mat=presence_matrix_filt,
                       #type="spp_site",
                       type="site_spp",
                       thresh=TRUE,
                       spp_names=TRUE)

cooccur_out <- cooccur_parallel(mat=tmp,
                                #type="spp_site",
                                type="site_spp",
                                thresh=TRUE,
                                spp_names=TRUE)

cooccur_parallel

summary(cooccur_out)

plot(cooccur_out)







presence_matrix_filt <- t(presence_matrix_filt)
# Note that this is JSD, not jaccard, so it's invalid!!
jaccard_net <- netConstruct(presence_matrix_filt,
                           filtTax = "none",
                           filtSamp = "none",
                           measure = "jsd",
                           normMethod = "none", 
                           zeroMethod = "none",
                           sparsMethod = "threshold",
                           thresh = 0.5,
                           verbose = 2,
                           seed = 131)

props_jaccard_net <- netAnalyze(jaccard_net, 
                           centrLCC = TRUE,
                           clustMethod = "cluster_fast_greedy",
                           hubPar = "eigenvector",
                           weightDeg = FALSE,
                           normDeg = FALSE)

p <- plot(props_jaccard_net, 
          nodeColor = "cluster", 
          nodeSize = "eigenvector",
          title1 = "Jaccard Network", 
          showTitle = TRUE,
          cexTitle = 2.3)

legend(0.7,
       1.1,
       cex = 2.2,
       title = "estimated association:",
       legend = c("+","-"),
       lty = 1,
       lwd = 3,
       col = c("#009900"
               "red"), 
       bty = "n",
       horiz = TRUE)


jaccard_net_gephi <- net_to_gephi(jaccard_net)

# Export to Gephi.
# (Taken from NetCoMi tutorial).
# Create edge object from the edge list exported by netConstruct()

net_to_gephi <- function(microNet_in, out_prefix = NULL) {
  
  edges <- dplyr::select(microNet_in$edgelist1, v1, v2)
  
  # Add Source and Target variables (as IDs)
  edges$Source <- as.numeric(factor(edges$v1))
  edges$Target <- as.numeric(factor(edges$v2))
  edges$Type <- "Undirected"
  edges$Weight <- microNet_in$edgelist1$adja
  
  nodes <- unique(edges[,c('v1','Source')])
  colnames(nodes) <- c("Label", "Id")
  
  # Add category with clusters (can be used as node colors in Gephi)
  nodes$Category <- microNet_in$clustering$clust1[nodes$Label]
  
  edges <- dplyr::select(edges, Source, Target, Type, Weight)
  
  if (! is.null(out_prefix)) {
    write.csv(nodes, file = paste(out_prefix, "_nodes.csv", sep = ""), row.names = FALSE)
    write.csv(edges, file = paste(out_prefix, "_edges.csv", sep = ""), row.names = FALSE)
  }
  
  return(list(edges=edges,
              nodes=nodes))
}
