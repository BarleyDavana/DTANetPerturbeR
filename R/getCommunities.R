#' Identify modules in a disease network
#'
#' @param disease_net A list of edges representing a disease network
#'
#' @return A list containing two data frames: "nodes" with two columns "name" and "module",
#' where "name" is the name of the node and "module" is the index of the module to which the node belongs;
#' and "edges" with three columns "name1", "name2", and "module", where "name1" and "name2" are the names
#' of the nodes connected by the edge and "module" is the index of the module to which the edge belongs.
#' Both data frames are saved as text files in the "output_Files" directory with the names "node_Module.txt"
#' and "edge_Module.txt", respectively.
#'
#' @export
#'
#' @examples
#' modules_info <- getCommunities(disease_net)

getCommunities <- function(disease_net) {
  library(igraph)
  set.seed(123)

  # Convert the disease network to a graph object
  edges <- data.frame(from = disease_net[, 1], to = disease_net[, 2])
  net <- graph_from_data_frame(edges, directed = FALSE)

  # Perform community detection using leading eigenvector and Louvain methods
  cle <- cluster_leading_eigen(net, weights = NULL)
  clv <- cluster_louvain(net, weights = NULL)

  # Associate each node with its corresponding community
  cle_label <- data.frame(label = get.vertex.attribute(net)[[1]],
                          community = cle$membership)
  clv_label <- data.frame(label = get.vertex.attribute(net)[[1]],
                          community = clv$membership)
  cle_label_list <- split(cle_label, f = cle_label$community)

  # Split nodes into separate lists based on their community assignment
  cle_label_list <- split(cle_label$label, f = cle_label$community)
  clv_label_list <- split(clv_label$label, f = clv_label$community)

  # Remove communities with only one node
  cle_label_list2 <- cle_label_list[lengths(cle_label_list) > 1]
  clv_label_list2 <- clv_label_list[lengths(clv_label_list) > 1]

  # Calculate the hypergeometric p-values between each pair of communities
  cle_LP2 <- as.data.frame(matrix(nrow = length(clv_label_list2), ncol = length(cle_label_list2)))
  cle_LP_union2 <- cle_LP_phyper2 <- cle_LP2
  for(i in 1:length(cle_label_list2)){
    for(j in 1:length(clv_label_list2)){
      cle_LP2[j,i]=length(intersect(cle_label_list2[[i]],
                                    clv_label_list2[[j]]))
      cle_LP_union2[j,i]=length(union(cle_label_list2[[i]],
                                      clv_label_list2[[j]]))
    }
  }
  for(i in 1:ncol(cle_LP_phyper2)){
    for(j in 1:nrow(cle_LP_phyper2)){
      cle_LP_phyper2[j,i]=1-phyper(cle_LP2[j,i],
                                   length(cle_label_list2[[i]]),cle_LP_union2[j,i],
                                   length(clv_label_list2[[j]]))
    }
  }

  # Assign module numbers to nodes based on significant p-values
  colnames(cle_LP_phyper2)<-paste("GN",
                                  seq(from=1,to=length(cle_label_list2),by=1),sep="_")
  rownames(cle_LP_phyper2)<-paste("LP",
                                  seq(from=1,to=length(clv_label_list2),by=1),sep="_")
  sig <- which(cle_LP_phyper2 < 0.05, arr.ind = TRUE)
  sig_list <- apply(sig, 1, function(x) intersect(cle_label_list2[[x[2]]], clv_label_list2[[x[1]]]))

  # Calculate the total number of nodes in the sig_list
  w <- sum(lengths(sig_list))
  modules <- data.frame(matrix(nc = 2, nr = w))

  # Create a new data frame with node and module columns
  modules <- do.call(rbind, Map(function(i, x) {
    data.frame(node = unlist(x), module = i)
  }, seq_along(sig_list), sig_list))

  # Subset the edges data frame to only include the first two columns and add a module column
  edges=edges[,1:2]
  edges$module <- 0

  # Iterate through each element in sig_list and each row in the edges data frame
  for(i in 1:length(sig_list))
    for(j in 1:nrow(edges))
    {
      pattern<-as.vector(modules[which(modules[,2]==i),1])
      if((edges[j,1] %in% pattern) &
         (edges[j,2] %in% pattern))
        edges$module[j] = i
    }

  # Subset the edges data frame to only include the rows with a module number assigned
  patternM<-edges[!edges$module==0,]
  colnames(modules) <- c("name","module")
  colnames(patternM) <- c("name1","name2","module")
  modules_info <- list(nodes = modules, edges = patternM)

  # Write the modules and edges data frames to separate output files
  write.table(modules, file = "output_Files/node_Module.txt",
              sep = "\t", col.names = T, row.names = F, quote = F)
  write.table(patternM, file = "output_Files/edge_Module.txt",
              sep = "\t", col.names = T, row.names = F, quote = F)

  return(modules_info)
}
