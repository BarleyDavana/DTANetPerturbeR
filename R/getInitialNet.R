#' Generate a sub-network consisting of the common genes and their interacting partners
#'
#' This function takes in a dataset of common genes and an output network name, and returns a sub-network that includes
#' only the common genes and their interacting partners from the human protein-protein interaction network. The sub-network
#' and a dataset of seed genes are also written to separate files in the output_Files folder.
#'
#' @param common_genes A dataset containing a list of common genes
#' @param outnet_name The name to be given to the output network files
#'
#' @return A list containing the dataset of seed genes and the sub-network
#' @export
#'
#' @examples
#' common_genes <- getCommonGenes(InterestedDGs)
#' seed_info <- getInitialNet(common_genes, "Initial_Network")
#' seed_genes <- seed_info$seed_genes
#' seed_net <- seed_info$seed_net
#' # The output files are saved in the output_Files folder as "Initial_Network_Genes.txt" and "Initial_Network_Net.txt"
#' # respectively.

getInitialNet <- function(common_genes, outnet_name){
  library(igraph)

  common_genes <- common_genes[, 1]
  data("PPIN_human")

  # Convert the edgelist to a graph object and get the nodes
  edges <- data.frame(from = PPIN_human[, 1], to = PPIN_human[, 2])
  net <- graph_from_data_frame(edges, directed = FALSE)
  nodes <- get.vertex.attribute(net)[[1]]

  # Get the sub-network containing only the common genes and their neighbors
  seeds <- intersect(nodes, common_genes)
  sub_net <- induced_subgraph(net, seeds)
  seed_nodes <- V(sub_net)[degree(sub_net) > 0]$name
  sub_net <- induced_subgraph(net, seed_nodes)

  # Get a dataset of seed genes and the sub-network as an edgelist
  seed_genes <- get.vertex.attribute(sub_net)[[1]]
  seed_genes <- data.frame(Gene_Symbol = seed_genes)
  seed_net <- get.edgelist(sub_net)
  colnames(seed_net) <- c("protein1", "protein2")
  seed_info <- list(seed_genes = seed_genes, seed_net = seed_net)

  # Write the seed genes and sub-network to separate files
  write.table(seed_genes, file = paste("output_Files/", outnet_name,
                                       "_nodes.txt", sep=""),
              sep = "\t", col.names = T, row.names = F, quote = F)
  write.table(seed_net, file = paste("output_Files/", outnet_name,
                                     "_edges.txt", sep=""),
              sep = "\t", col.names = T, row.names = F, quote = F)

  return(seed_info)
}
