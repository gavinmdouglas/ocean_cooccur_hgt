rm(list = ls(all.names = TRUE))

tmp <- readRDS("~/Desktop/working/set1.rds")

library(igraph)


tmp_graph <- igraph::graph_from_data_frame(d = tmp[, c("taxoni", "taxonj", "simple_cooccur")],
                                           directed = FALSE, vertices = NULL)

igraph::graph_from_edgelist(directed = FALSE)

igraph::plot.igraph()

plot(tmp_graph, vertex.label=NA, vertex.size=0.5)


igraph::write_graph(tmp_graph, file="~/Desktop/working/set1.graphml2",
                    format="graphml")


degree_breakdown <- igraph::degree(tmp_graph)
closeness_centrality <- igraph::closeness(tmp_graph, mode="all")
betweenness_centrality <- igraph::betweenness(tmp_graph, directed = FALSE)