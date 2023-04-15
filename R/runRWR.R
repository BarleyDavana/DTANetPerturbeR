#' Run Random Walk with Restart (RWR) on a protein-protein interaction network (PPIN)
#'
#' This function runs RWR on a PPIN to prioritize genes based on their network proximity to the input seeds.
#'
#' @param seeds A data frame containing a list of gene symbols to be used as input seeds for RWR.
#' @param restart A numeric value indicating the restart probability for RWR (default = 0.8).
#' @param output A logical value indicating whether or not to output the RWR results as a text file (default = FALSE).
#'
#' @return A data frame containing the RWR scores for each gene in the PPIN.
#' @export
#'
#' @examples
#' # Run RWR with default restart probability and no output file
#' runRWR(seeds = data.frame("BRCA1", "TP53"))
#'
#' # Run RWR with custom restart probability and output file
#' runRWR(seeds = data.frame("BRCA1", "TP53"), restart = 0.7, output = TRUE)

runRWR <- function(seeds, restart = 0.8, output = FALSE){
  library(igraph)
  library(dnet)
  library(dplyr)

  # Load human protein-protein interaction network data
  data("PPIN_human")

  # Convert edge data frame to an undirected graph object
  edges <- data.frame(from = PPIN_human[, 1], to = PPIN_human[, 2])
  net <- graph.data.frame(edges, directed = FALSE)
  seeds <- seeds[, 1]

  # Assign labels to the nodes in the network
  nodes <- get.vertex.attribute(net)[[1]]
  nodes <- as.data.frame(nodes)
  rownames(nodes) <- nodes$nodes
  nodes$selected <- 0
  nodes[seeds, 2] <- 1
  nodes <- nodes[order(nodes[, 2], decreasing = T), ]

  # Count the number of seed nodes
  seedcount <- as.numeric(length(which(nodes[, 2] == 1)))
  g <- graph_from_data_frame(d = edges, vertices = nodes, directed = F)

  # Initialize a matrix to store the RWR scores for each node
  Seeds <- as.data.frame(matrix(nc = 1, nr = vcount(net)))
  Seeds$V1 <- 0
  rownames(Seeds) <- nodes[, 1]
  Seeds[1:seedcount, 1] <- 1
  result <- as.data.frame(matrix(nrow = vcount(net), ncol = seedcount))

  # Compute RWR scores for each seed node
  for(i in 1:seedcount){
    Seeds[1:seedcount, 1] <- 0
    Seeds[i, 1] <- 1
    PTmatrix <- dRWR(g = g, normalise = "laplacian",
                     setSeeds = Seeds, restart = restart, parallel = FALSE)
    result[, i] <- PTmatrix[, 1]
  }
  colnames(result) <- nodes[1:seedcount, 1]
  rownames(result) <- nodes[, 1]

  # Compute total RWR score for each node if there are multiple seed nodes
  if (seedcount != 1) {
    result$RWRs <- apply(result[, c(rep(1:seedcount, 1))], 1, sum, na.rm=T)
  } else {
    result$RWRs <- result[, 1]
  }
  result <- mutate(result, rn = row.names(result))

  # Create a data frame of RWR scores for each node
  RWR_score <- as.data.frame(result$rn)
  RWR_score$RWRs <- result$RWRs
  RWR_score <- arrange(RWR_score, desc(RWRs))
  colnames(RWR_score) <- c("Gene_Symbol", "RWR_score")
  if (output){
    write.table(RWR_score, paste("output_Files/", 'RWR_',
                                 restart, '_score.txt', sep = ''),
                quote = F, sep = "\t", row.names = F)
  }

  return(RWR_score)
}
