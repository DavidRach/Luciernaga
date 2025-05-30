% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/QC_SimilarFluorophores.R
\name{QC_SimilarFluorophores}
\alias{QC_SimilarFluorophores}
\title{Queries fluorophore and returns similar fluorophores.}
\usage{
QC_SimilarFluorophores(
  TheFluorophore,
  NumberDetectors,
  returnSynonyms = FALSE,
  NumberHits = 10,
  returnPlots = FALSE,
  plotlinecolor = NULL
)
}
\arguments{
\item{TheFluorophore}{The name of the Fluorophore compare, see QC_ReferenceLibrary}

\item{NumberDetectors}{Number of detectors of the instrument}

\item{returnSynonyms}{Something}

\item{NumberHits}{Number of most similar fluorophores by cosine}

\item{returnPlots}{Whether to also return signature plots, default is set FALSE}

\item{plotlinecolor}{Default NULL, otherwise if single line provide desired color}

\item{returnSynonymns}{Returns only fluorophores > 0.98 cosine value, default FALSE}
}
\value{
A dataframe of similar fluorophores
}
\description{
Queries fluorophore and returns similar fluorophores.
}
\examples{
Results <- QC_SimilarFluorophores(TheFluorophore="Spark Blue 550",
 NumberDetectors=64, returnSynonymns=FALSE, NumberHits = 10, returnPlots=FALSE)
}
