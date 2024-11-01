% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/QC_ReferenceLibrary.R
\name{QC_ReferenceLibrary}
\alias{QC_ReferenceLibrary}
\title{Queries the available reference library for available fluorophores (and their naming conventions)}
\usage{
QC_ReferenceLibrary(FluorNameContains, NumberDetectors, returnPlots = FALSE)
}
\arguments{
\item{FluorNameContains}{A character string pattern to match, example "APC"}

\item{NumberDetectors}{The Number of Detectors for your instrument}

\item{returnPlots}{Whether to return signature plot as well. Default FALSE.}
}
\value{
A dataframe column containing matching Fluorophores from your querry
}
\description{
Queries the available reference library for available fluorophores (and their naming conventions)
}
\examples{
QC_ReferenceLibrary(FluorNameContains = "FITC", NumberDetectors=64)
}
