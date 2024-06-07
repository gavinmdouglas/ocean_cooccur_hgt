rm(list = ls(all.names = TRUE))

library(NetCoMi)
library(SpiecEasi)

metaG_rpkm <- read.table("/mfs/gdouglas/projects/ocean_mags/coverm/combined_tables/metaG_rpkm_alldata.tsv.gz",
                         header = TRUE, sep = "\t", row.names = 1)

metaG_rpkm_netcomi <- netConstruct(metaG_rpkm,
                                   measure = "spieceasi",
                                   verbose = 3,
                                   seed = 123456)

colnames(metaG_rpkm_netcomi$edgelist1)[1] <- "taxon_i"
colnames(metaG_rpkm_netcomi$edgelist1)[2] <- "taxon_j"

write.table(x = metaG_rpkm_netcomi$edgelist1,
            file = "/mfs/gdouglas/projects/ocean_mags/coverm/network_working_alldata/metaG_spieceasi_rpkm.alldata.tsv",
            sep = "\t",
            quote = FALSE,
            col.names = TRUE,
            row.names = FALSE)
