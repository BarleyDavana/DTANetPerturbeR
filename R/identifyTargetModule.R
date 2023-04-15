#' Identify target module in disease gene network
#'
#' This function first detects modules in a given disease gene network, then maps targets to genes
#' in the network. It then identifies the target module, and maps drugs to targets in the target module.
#' The final result is a data frame of drug-target interactions in the target module.
#'
#' @param disease_net A list of edges representing a disease network
#'
#' @return NULL
#' @export
#'
#' @examples
#' # Creating a disease network using createDiseaseNetwork function
#' disease_net <- createDiseaseNetwork(gene_list, enlarge_network = FALSE)
#'
#' # Creating a disease network from a text file containing network edges
#' disease_net <- read.table("disease_net.txt", header = TRUE)
#'
#' # Identifying the target module in the network using identifyTargetModule function
#' identifyTargetModule(disease_net)

identifyTargetModule <- function(disease_net) {
  # Create output directory
  if (!file.exists("output_Files")) {
    dir.create("output_Files")
  }

  # Step 1: Module detection
  modules_info <- getCommunities(disease_net)

  # Step 2: Map targets to genes in disease net
  net_tars <- mapTargets(disease_net)

  # Step 3: Identify target module
  target_module_info <- getTargetModule(modules_info, net_tars)

  # Step 4: Map drugs to targets in target module
  mapDrugs(target_module_info)
}
