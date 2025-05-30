% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/QC_GainMonitoring.R
\name{QC_GainMonitoring}
\alias{QC_GainMonitoring}
\title{Takes Before After QC .fcs files of QC beads run on assay settings, returns the gains and average
measurement in a data.frame to allow for plotting to assess stability over time.}
\usage{
QC_GainMonitoring(
  x,
  sample.name,
  stats,
  subsets = NULL,
  inverse.transform = FALSE
)
}
\arguments{
\item{x}{A CytoSet (when mapped) or an individual cytoset, example (x=MyCytoSet[[1]])}

\item{sample.name}{The keyword value that distinguishes individual .fcs files}

\item{stats}{Whether to summarize the MFIs by "mean" or "median"}

\item{subsets}{When provided with a GatingSet, selects this subset to provide the Cytoset}

\item{inverse.transform}{Default is FALSE.}
}
\value{
A data.frame row
}
\description{
Takes Before After QC .fcs files of QC beads run on assay settings, returns the gains and average
measurement in a data.frame to allow for plotting to assess stability over time.
}
\examples{
library(dplyr)
library(purrr)
library(flowWorkspace)

File_Location <- system.file("extdata", package = "Luciernaga")
FCS_Pattern <- ".fcs$"
FCS_Files <- list.files(path = File_Location, pattern = FCS_Pattern,
full.names = TRUE, recursive = FALSE)
QCBeads <- FCS_Files[grep("Before|After", FCS_Files)]
BeforeAfter_CS <- load_cytoset_from_fcs(files=QCBeads,
transformation=FALSE, truncate_max_range = FALSE)
BeforeAfter <- map(.x=BeforeAfter_CS[1:2], .f=QC_GainMonitoring,
 sample.name = "TUBENAME", stats="median") \%>\% bind_rows()
}
