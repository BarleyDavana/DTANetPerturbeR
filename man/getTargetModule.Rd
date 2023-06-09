% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getTargetModule.R
\name{getTargetModule}
\alias{getTargetModule}
\title{Get information of the module with the highest proportion of target genes}
\usage{
getTargetModule(modules_info, Net_Tars)
}
\arguments{
\item{modules_info}{A list containing information on the modules, including node and edge information.}

\item{Net_Tars}{A data frame containing all the target genes in disease net.}
}
\value{
A list containing information on the target module, including its nodes, edges, and targets.
}
\description{
This function creates a module network graph for each module in the input data,
computes the proportion of target genes in each module, and selects the module with
the highest proportion of target genes. It then retrieves information on the selected
module, including its nodes, edges, and targets, and writes the edge information to a file.
}
\examples{
getTargetModule(modules_info, Net_Tars)
}
