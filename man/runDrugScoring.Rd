% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/runDrugScoring.R
\name{runDrugScoring}
\alias{runDrugScoring}
\title{runDrugScoring Function}
\usage{
runDrugScoring(python_exe, dti_file, fasta_file, network_file, output_folder)
}
\arguments{
\item{python_exe}{The path to the Python executable file.}

\item{dti_file}{The path to the input drug-target interaction (DTI) file in TSV format. The file should contain two columns: drug IDs and protein UniProt IDs.}

\item{fasta_file}{The path to the input protein FASTA file.}

\item{network_file}{The path to the input protein-protein interaction network file.}

\item{output_folder}{The path to the output directory where the results will be stored.}
}
\value{
This function doesn't return anything, but saves the results to the specified output folder.
}
\description{
This function runs two Python scripts for predicting drug-target binding affinity
and calculating a protein ranking score based on the predicted binding affinities.
}
\examples{
runDrugScoring(python_exe = "python", dti_file = "Target_module_dti.txt", fasta_file = "proteins.fasta",
network_file = "Target_module_edges.txt", output_folder = "output_Files")
}
