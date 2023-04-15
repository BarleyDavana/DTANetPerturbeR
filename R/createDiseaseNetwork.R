#' Create Disease Network
#'
#' This function creates a disease gene network for a given set of disease genes.
#' It first obtains the common genes shared among other diseases, then uses these
#' genes to construct an initial network. If the "enlarge_network" parameter is set
#' to TRUE, the function enlarges the network. The final disease network is returned.
#'
#' @param disease_genes A data frame containing disease gene names in a single column.
#' @param enlarge_network A logical parameter to indicate whether to enlarge the network or not.
#'
#' @return The final disease gene network.
#' @export
#'
#' @examples
#' # Load interested disease genes from file
#' gene_list <- read.table("gene_list.txt", header = T, sep = '\t')
#'
#' # Create the disease network
#' createDiseaseNetwork(gene_list)

createDiseaseNetwork <- function(disease_genes, enlarge_network = FALSE){
  # Create output directory
  if (!file.exists("output_Files")) {
    dir.create("output_Files")
  }

  # Step 1: Get common genes with other diseases
  common_genes <- getCommonGenes(disease_genes)

  # Step 2: Get initial network using common genes
  seed_info <- getInitialNet(common_genes, "Initial_Network")
  seed_genes <- seed_info$seed_genes
  seed_net <- seed_info$seed_net

  # Step 3: Enlarge network if enlarge_network parameter is TRUE
  if (enlarge_network){
    disease_net <- getEnlargedNet(seed_genes)
  } else {
    disease_net <- seed_net
  }

  # Step 4: Output final disease network
  return(disease_net)
}
