rm(list = ls(all.names = TRUE))

library(NetCoMi)

metaG_rpkm <- read.table("/mfs/gdouglas/projects/water_mags/coverm/combined_tables/metaG_rpkm.tsv.gz",
                         header = TRUE, sep = "\t", row.names = 1)

metaT_rpkm <- read.table("/mfs/gdouglas/projects/water_mags/coverm/combined_tables/metaT_rpkm.tsv.gz",
                         header = TRUE, sep = "\t", row.names = 1)

metaG_rpkm_netcomi <- netConstruct(metaG_rpkm,
                                   measure = "propr",
                                   sparsMethod="threshold",
                                   thresh = 0.3,
                                   zeroMethod="multRepl",
                                   verbose = 3,
                                   seed = 123456)

colnames(metaG_rpkm_netcomi$edgelist1)[1] <- "taxon_i"
colnames(metaG_rpkm_netcomi$edgelist1)[2] <- "taxon_j"

write.table(x = metaG_rpkm_netcomi$edgelist1,
            file = "/mfs/gdouglas/projects/water_mags/coverm/network_working/metaG_propr_rpkm.tsv",
            sep = "\t",
            quote = FALSE,
            col.names = TRUE,
            row.names = FALSE)


metaT_rpkm_netcomi <- netConstruct(metaT_rpkm,
                                   measure = "propr",
                                   sparsMethod="threshold",
                                   thresh = 0.3,
                                   zeroMethod="multRepl",
                                   verbose = 3,
                                   seed = 123456)

colnames(metaT_rpkm_netcomi$edgelist1)[1] <- "taxon_i"
colnames(metaT_rpkm_netcomi$edgelist1)[2] <- "taxon_j"

write.table(x = metaT_rpkm_netcomi$edgelist1,
            file = "/mfs/gdouglas/projects/water_mags/coverm/network_working/metaT_propr_rpkm.tsv",
            sep = "\t",
            quote = FALSE,
            col.names = TRUE,
            row.names = FALSE)
