#' runDrugScoring Function
#'
#' This function runs two Python scripts for predicting drug-target binding affinity
#' and calculating a protein ranking score based on the predicted binding affinities.
#'
#' @param python_exe The path to the Python executable file.
#' @param dti_file The path to the input drug-target interaction (DTI) file in TSV format. The file should contain two columns: drug IDs and protein UniProt IDs.
#' @param fasta_file The path to the input protein FASTA file.
#' @param network_file The path to the input protein-protein interaction network file.
#' @param output_folder The path to the output directory where the results will be stored.
#'
#' @return This function doesn't return anything, but saves the results to the specified output folder.
#'
#' @export
#'
#' @examples
#' runDrugScoring(python_exe = "python", dti_file = "Target_module_dti.txt", fasta_file = "proteins.fasta",
#' network_file = "Target_module_edges.txt", output_folder = "output_Files")

runDrugScoring <- function(python_exe, dti_file, fasta_file, network_file, output_folder) {
  # Create output directory
  if (!file.exists("output_Files")) {
    dir.create("output_Files")
  }

  # Step 1: Set paths to Python scripts
  predict_binding_affinity <- system.file("predict_binding_affinity.py", package = "DTANetPerturbeR")
  calculate_PRS <- system.file("calculate_PRS.py", package = "DTANetPerturbeR")

  # Step 2: Run dti_prediction.py
  system2(command = python_exe, args = c(predict_binding_affinity, shQuote(dti_file), shQuote(output_folder), shQuote(fasta_file)))

  # Step 3: Run calculate_PRS.py
  affinity_file <- file.path(output_folder, "virtual_screening.txt")
  system2(command = python_exe, args = c(calculate_PRS, shQuote(network_file), shQuote(affinity_file), shQuote(output_folder)))
}
