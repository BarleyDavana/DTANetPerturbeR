#' Map drugs to target module
#'
#' @param target_module_info list of nodes, edges, and targets for the target module
#'
#' @return NULL
#'
#' @export
#'
#' @examples
#' mapDrugs(target_module_info)

mapDrugs <- function(target_module_info){
  library(data.table)
  library(dplyr)
  library(clusterProfiler)
  library(org.Hs.eg.db)

  data("Drugs_DrugBank")
  data("DTI_DrugBank")
  data("Drugs_TTD")
  data("Drugs_supp_TTD")
  data("DTI_TTD")

  # Extract target information from input
  targets <- target_module_info$targets
  t <- bitr(targets, fromType="SYMBOL", toType="UNIPROT", OrgDb="org.Hs.eg.db")

  # Subset DrugBank DTI data for UniProt IDs present in target list
  dti <- DTI_DrugBank[DTI_DrugBank$`UniProt ID` %in% t$UNIPROT, ]
  result <- inner_join(t, dti, by = c("UNIPROT" = "UniProt ID"))
  dti_DrugBank <- inner_join(result, Drugs_DrugBank)
  dti_DrugBank <- dti_DrugBank[, c(1:4, ncol(dti_DrugBank))]
  colnames(dti_DrugBank) <- c("GENENAME", "UNIPROTID", "DRUGID", "DRUGNAME", "SMILES")
  dti_DrugBank$DATABASE <- rep("DrugBank", nrow(dti_DrugBank))

  # Extract drug SMILES from TTD data
  DRUGSMIL <- Drugs_TTD %>%
    dplyr::filter(V2 == "DRUGSMIL") %>%
    dplyr::select(V1, V3) %>%
    dplyr::rename(DrugID = V1, DRUGSMIL = V3)

  # Extract drug names from TTD supplemental data
  DRUGNAME <- Drugs_supp_TTD %>%
    dplyr::filter(V2 == "DRUGNAME") %>%
    dplyr::select(V1, V3) %>%
    dplyr::rename(DrugID = V1, DRUGNAME = V3)

  # Merge drug SMILES and names from TTD data
  drugs_TTD <- inner_join(DRUGSMIL, DRUGNAME)

  # Extract target gene names from TTD data
  targets_TTD <- Tars_TTD %>%
    dplyr::filter(V2=="GENENAME") %>%
    dplyr::select(V1, V3) %>%
    dplyr::rename(TargetID = V1, GENENAME = V3)

  # Merge drug, target, and interaction data from TTD
  dti <- DTI_TTD %>%
    inner_join(drugs_TTD) %>%
    inner_join(targets_TTD)
  t <- bitr(targets, fromType="SYMBOL", toType="UNIPROT", OrgDb="org.Hs.eg.db")
  colnames(t) <- c("GENENAME", "UNIPROT")
  dti_TTD <- inner_join(t, dti)

  # Subset the columns
  dti_TTD <- dti_TTD[, c(1, 2, 4, 7, 8)]
  colnames(dti_TTD) <- c("GENENAME", "UNIPROTID", "DRUGID", "SMILES", "DRUGNAME")
  dti_TTD$DATABASE <- rep("TTD", nrow(dti_TTD))

  # Combine the DTI data from both DrugBank and TTD
  Target_module_dti <- full_join(dti_DrugBank, dti_TTD)

  # Write the final merged dataframe to a text file
  write.table(Target_module_dti, file = "output_Files/Target_module_dti.txt",
              sep = "\t", col.names = T, row.names = F, quote = F)
}
