% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/createDiseaseNetwork.R
\name{createDiseaseNetwork}
\alias{createDiseaseNetwork}
\title{Create Disease Network}
\usage{
createDiseaseNetwork(disease_genes, enlarge_network = FALSE)
}
\arguments{
\item{disease_genes}{A data frame containing disease gene names in a single column.}

\item{enlarge_network}{A logical parameter to indicate whether to enlarge the network or not.}
}
\value{
The final disease gene network.
}
\description{
This function creates a disease gene network for a given set of disease genes.
It first obtains the common genes shared among other diseases, then uses these
genes to construct an initial network. If the "enlarge_network" parameter is set
to TRUE, the function enlarges the network. The final disease network is returned.
}
\examples{
# Load interested disease genes from file
gene_list <- read.table("gene_list.txt", header = T, sep = '\t')

# Create the disease network
createDiseaseNetwork(gene_list)
}
