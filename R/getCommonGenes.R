#' Get common genes
#'
#' This function intersects a dataset of unique gene symbols with a dataset of interested genes,
#' and returns only the genes that are common to both datasets,
#' representing the genes that are shared between the disease of interest and other diseases.
#'
#' @param InterestedDGs A dataset containing a list of genes of interest
#'
#' @return A dataset with only the common genes
#' @export
#'
#' @examples
#' # Load interested disease genes from file
#' gene_list <- read.table("gene_list.txt", header = T, sep = '\t')
#'
#' # Get common genes with a gene dataset
#' getCommonGenes(gene_list)

getCommonGenes <- function(InterestedDGs){
  library(limma)

  # Process the disease genes dataset to get a dataset of unique gene symbols
  data("Disease_Genes")
  DGs <- Disease_Genes[, 1]

  # Split disease genes by comma
  SplitDGs <- c()
  for(i in 1:length(DGs)){
    split <- strsplit2(DGs[i], ",")
    split2 <- as.vector(split)
    SplitDGs <- c(SplitDGs, split2)
  }
  UniqueDGs <- SplitDGs[!duplicated(SplitDGs)]

  # Find the common genes between UniqueDGs and InterestedDGs
  InterestedDGs <- InterestedDGs[, 1]
  common_genes <- intersect(UniqueDGs, InterestedDGs)
  common_genes <- data.frame(common_genes = common_genes)

  # Write the common genes to a file
  write.table(common_genes, file = "output_Files/Common_Genes.txt",
              sep = "\t", col.names = T, row.names = F, quote = F)

  return(common_genes)
}
