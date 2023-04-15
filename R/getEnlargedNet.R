#' Extracts an enlarged subgraph of a given seed network based on reachability score
#'
#' This function takes a seed network, calculates reachability scores of each non-seed gene using
#' Random Walk with Restart (RWR) algorithm, calculates fold change and p-values between seed and non-seed genes,
#' identifies significant genes based on the p-value and log2 fold change values, and then extracts an enlarged
#' subgraph that includes the seed genes and the significant non-seed genes.
#' The enlargement threshold can be adjusted using the reachable_genes_threshold parameter.
#'
#' @param seed_genes a data frame containing the seed genes
#' @param reachable_genes_threshold a numeric value representing the reachability score threshold for identifying significant genes
#'
#' @return a data frame containing the enlarged subgraph of the seed network
#' @export
#'
#' @examples
#' enlarged_net <- getEnlargedNet(seed_genes, 0.9)

getEnlargedNet <- function(seed_genes, reachable_genes_threshold = 0.9){
  # Run RWR algorithm on the seed genes
  seedRWR <- runRWR(seed_genes)

  # Get the list of seed genes
  seedgenes <- seed_genes[, 1]
  seedcount <- length(seedgenes)

  # Separate the nodes into seed and non-seed nodes
  nodes <- seedRWR[, 1]
  seeds <- nodes[1:seedcount]
  non_seeds <- nodes[(seedcount+1):length(nodes)]

  # Initialize empty vectors for Fold Change and p-value
  Fold_Change <- c()
  Pvalue <- c()

  # Loop through each non-seed node
  for (non_seed in non_seeds) {
    node <- data.frame(node = non_seed)

    # Run RWR algorithm on the non-seed node
    reachability_score <- runRWR(node)
    rs <- reachability_score[-1, ]

    # Get the list of reachable genes from seed nodes and from non-seed nodes
    reach_list_seeds <- rs[rs[, 1] %in% seeds, 2]
    reach_list_non_seeds <- rs[rs[, 1] %in% non_seeds, 2]
    pvalue <- wilcox.test(reach_list_seeds, reach_list_non_seeds,
                          alternative = "two.sided")$p.value
    Pvalue <- c(Pvalue, pvalue)

    # Calculate the Fold Change score
    reachable_seeds <- mean(rs[rs[, 1] %in% seeds, 2])
    reachable_non_seeds <- mean(rs[rs[, 1] %in% non_seeds, 2])
    FC <- reachable_seeds/reachable_non_seeds
    Fold_Change <- c(Fold_Change, FC)
  }

  # Calculate the log2 Fold Change and store the results in a data frame
  Fold_Change_log <- log2(Fold_Change)
  result <- data.frame(non_seeds, Fold_Change,
                       Fold_Change_log, Pvalue)
  colnames(result) <- c("Gene_Symbol", "Fold_Change", "log(FC)", "Pvalue")

  # Write the result to a file
  write.table(result, "output_Files/Reachability_Score.txt",
              quote = F, sep = "\t", row.names = F)
  Reachability <- result

  # Filter the significant genes based on p-value and log2 Fold Change
  sig_genes <- Reachability$Pvalue < 0.05
  largefc_genes <- Reachability$`log(FC)` > 2
  reachable_genes <- sig_genes & largefc_genes
  selected_genes <- Reachability$Gene_Symbol[reachable_genes]
  rwr_rank <- which(Reachability$Gene_Symbol %in% selected_genes)
  reachable_genes <- data.frame(rwr_rank, selected_genes)

  # Mark the selected genes in the Reachability data frame
  Reachability$ifselected <- 0
  selected_genes <- reachable_genes$selected_genes
  rownames(Reachability) <- Reachability$Gene_Symbol
  Reachability[selected_genes, 5] <- 1

  # Calculate the cumulative number of selected genes
  Reachability$genenum <- 1:nrow(Reachability)
  Reachability$selected_genes_num <- cumsum(Reachability$ifselected)
  Reachability$density <- Reachability$selected_genes_num/length(selected_genes)

  # Store the density values of reachable genes
  reachable_genes_dens <- Reachability$density

  # Find the index of the first row in Reachability matrix where the density is above the threshold
  sub_nodes_idx <- which.max(reachable_genes_dens >= reachable_genes_threshold)
  netsize <- sub_nodes_idx + 80

  # Select the subgraph up to the row determined by netsize
  sub_nodes <- seedRWR[1:netsize, ]
  Dnet_info <- getInitialNet(sub_nodes, "Enlarged_Network")
  Enlarged_net <- Dnet_info$seed_net

  return(Enlarged_net)
}
