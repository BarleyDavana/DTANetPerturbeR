#' Get information of the module with the highest proportion of target genes
#'
#' This function creates a module network graph for each module in the input data,
#' computes the proportion of target genes in each module, and selects the module with
#' the highest proportion of target genes. It then retrieves information on the selected
#' module, including its nodes, edges, and targets, and writes the edge information to a file.
#'
#' @param modules_info A list containing information on the modules, including node and edge information.
#' @param Net_Tars A data frame containing all the target genes in disease net.
#'
#' @return A list containing information on the target module, including its nodes, edges, and targets.
#' @export
#'
#' @examples
#' getTargetModule(modules_info, Net_Tars)

getTargetModule <- function(modules_info, Net_Tars){
  library(igraph)
  edges <- modules_info$edges

  # Get all of the modules
  all_modules <- unique(edges$module)

  # Initialize the maximum percentage of targets to 0
  max_perc_targets <- 0

  # Loop through each module
  for (m in all_modules) {
    module_edges <- edges[edges$module == m, ]
    module_edges <- data.frame(from = module_edges[, 1], to = module_edges[, 2])
    module_net <- graph_from_data_frame(module_edges, directed = FALSE)
    module_nodes <- get.vertex.attribute(module_net)[[1]]
    all_tars <- Net_Tars[, 1]
    module_tars <- intersect(all_tars, module_nodes)

    # Calculate the percentage of targets
    perc_targets <- length(module_tars)/length(module_nodes)
    # Update the target module information
    if (perc_targets > max_perc_targets) {
      max_perc_targets <- perc_targets
      target_module_edges <- module_edges
      target_module_nodes <- module_nodes
      target_module_targets <- module_tars
    }
  }

  # Create a network object for the target module
  target_module_net <- graph_from_data_frame(
    data.frame(from = target_module_edges[, 1],
               to = target_module_edges[, 2]),
    directed = FALSE)

  # Save the target module information in a list
  target_module_info <- list(nodes = target_module_nodes,
                             edges = target_module_edges,
                             targets = target_module_targets)

  # Save the target module edges in a file
  target_module_ppi <- target_module_edges[, 1:2]
  colnames(target_module_ppi) <- c('gene1', 'gene2')
  write.table(target_module_ppi, file = "output_Files/Target_module_edges.txt",
              sep = "\t", col.names = T, row.names = F, quote = F)

  return(target_module_info)
}
