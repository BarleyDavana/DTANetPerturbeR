#' Identify targets of a gene network using DrugBank and TTD databases
#'
#' @param Gene_net a data frame representing a gene network, where the first two columns represent
#' the gene symbols of the source and target nodes, respectively.
#'
#' @return a data frame containing the symbols of the genes in the input network that are also targets in
#' either DrugBank or TTD databases
#' @export
#'
#' @examples
#' targets <- mapTargets(disease_net)
#' head(targets)

mapTargets <- function(Gene_net){
  library(dplyr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(igraph)

  data("DTI_DrugBank")
  data("Tars_TTD")
  tars_drugbank <- DTI_DrugBank$`UniProt ID`
  tars_drugbank <- tars_drugbank[!duplicated(tars_drugbank)]

  # Convert UniProt IDs to gene symbols
  t = bitr(tars_drugbank, fromType="UNIPROT", toType="SYMBOL", OrgDb="org.Hs.eg.db")
  tars_drugbank <- t$SYMBOL[!duplicated(t$SYMBOL)]

  # Extract target genes from TTD database
  tars_TTD <- Tars_TTD %>%
    dplyr::filter(V2 == "GENENAME") %>%
    dplyr::select(V3) %>%
    dplyr::rename(SYMBOL = V3)
  tars_TTD <- tars_TTD$SYMBOL[!duplicated(tars_TTD$SYMBOL)]

  # Create a graph object from the input gene network
  edges <- data.frame(from = Gene_net[, 1], to = Gene_net[, 2])
  net <- graph.data.frame(edges, directed = FALSE)
  nodes <- get.vertex.attribute(net)[[1]]

  # Get target genes from DrugBank database that are also present in the input gene network
  net_tar_drugbank <- intersect(tars_drugbank, nodes)

  # Get target genes from TTD database that are also present in the input gene network
  net_tar_TTD <- intersect(tars_TTD, nodes)

  # Combine the target genes from both databases and write to file
  net_tars <- union(net_tar_drugbank, net_tar_TTD)
  net_tars <- data.frame(Gene_Symbol = net_tars)
  write.table(net_tars, file = "output_Files/Net_Tars.txt",
              sep = "\t", col.names = T, row.names = F, quote = F)

  return(net_tars)
}
