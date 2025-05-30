% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/QC_ViewSignature.R
\name{QC_ViewSignature}
\alias{QC_ViewSignature}
\title{Visualizes the Signature for given row in an averaged signature
 data.frame.}
\usage{
QC_ViewSignature(
  x,
  columnname = "Sample",
  data,
  Normalize = TRUE,
  TheFormat = "wider",
  detectorcolumn = NULL,
  valuecolumn = NULL,
  legend = TRUE,
  plotname = NULL,
  plotlinecolor = NULL
)
}
\arguments{
\item{x}{Name in the Sample column you want to filter for}

\item{columnname}{Default is Sample, specifies column name
from which x is filtered from}

\item{data}{A data.frame object from QC_LibraryParse containing
 Fluorophore name column 
and numeric detector columns.}

\item{Normalize}{Whether to normalize the data based on peak
detector value, default is TRUE}

\item{TheFormat}{Default wider for detectors in columns, specify
longer if providing detectors as rows}

\item{detectorcolumn}{Default NULL, when TheFormat="longer" specify
detector column name}

\item{valuecolumn}{Default NULL, when TheFormat="longer" specify
value column name}

\item{legend}{Default TRUE, alternately removes plot legend}

\item{plotname}{Default NULL, alternately specifies the plot
title}

\item{plotlinecolor}{Default NULL, alternatively provide color
when only a single line}
}
\value{
A dataframe of similar fluorophores
}
\description{
Visualizes the Signature for given row in an averaged signature
 data.frame.
}
\examples{

library(flowCore)
library(flowWorkspace)
library(openCyto)
library(data.table)
library(dplyr)

File_Location <- system.file("extdata", package = "Luciernaga")
FCS_Files <- list.files(path = File_Location, pattern = ".fcs",
  full.names = TRUE)
UnstainedFCSFiles <- FCS_Files[grep("Unstained", FCS_Files)]
UnstainedCells <- UnstainedFCSFiles[-grep("Beads", UnstainedFCSFiles)]
MyCytoSet <- load_cytoset_from_fcs(UnstainedCells[1],
  truncate_max_range = FALSE,transformation = FALSE)
MyGatingSet <- GatingSet(MyCytoSet)
MyGates <- fread(file.path(path = File_Location, pattern = 'Gates.csv'))
MyGatingTemplate <- gatingTemplate(MyGates)
gt_gating(MyGatingTemplate, MyGatingSet)
removestrings <-  c("DR_", "Cells", ".fcs", "-", " ")
StorageLocation <- file.path("C:", "Users", "JohnDoe", "Desktop")

PopulationInterest <- gs_pop_get_data(MyGatingSet[1], subset="lymphocytes")
TheDataValues <- exprs(PopulationInterest[[1]])
TheDataValues <- data.frame(TheDataValues, check.names=FALSE)
Signature <- AveragedSignature(TheDataValues, stats="median")
TheData <- Signature[,-grep("Time|FS|SC|SS|Original|W$|H$", names(Signature))]
TheData <- TheData |> mutate(Sample="TestSignature") |>
 relocate(Sample, .before=1)

Plot <- QC_ViewSignature(x="TestSignature", data=TheData, Normalize=TRUE)

}
