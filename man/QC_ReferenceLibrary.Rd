% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/QC_ReferenceLibrary.R
\name{QC_ReferenceLibrary}
\alias{QC_ReferenceLibrary}
\title{Queries the available reference library for available fluorophores (and their naming conventions)}
\usage{
QC_ReferenceLibrary(
  FluorNameContains,
  NumberDetectors,
  returnPlots = FALSE,
  ListOverride = FALSE
)
}
\arguments{
\item{FluorNameContains}{A character string pattern to match, example "APC"}

\item{NumberDetectors}{The number of detectors the instrument contains.
For Cytek Instruments 5L = 64, 4L_UV = 54, 4L_YG = 48, 3L=38, 2L_VB=30,
2L_BR=22, 1L=14}

\item{returnPlots}{Whether to return signature plot as well. Default FALSE.}

\item{ListOverride}{Default FALSE}
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
