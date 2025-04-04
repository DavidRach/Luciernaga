% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Utility_IterativeGating.R
\name{Utility_IterativeGating}
\alias{Utility_IterativeGating}
\title{Return individual plots for specified node for specimens within a GatingSet.}
\usage{
Utility_IterativeGating(
  x,
  subset,
  gate,
  xValue,
  yValue,
  sample.name,
  removestrings,
  bins = 270
)
}
\arguments{
\item{x}{A Gating Set Object}

\item{subset}{Desired node of cells to plot}

\item{gate}{Desired gate to show on the same plot as the subset, else NULL}

\item{xValue}{Desired x axis}

\item{yValue}{Desired y axis}

\item{sample.name}{Keyword under which the sample name is stored}

\item{removestrings}{Character values to remove from the name}

\item{bins}{Geo_hex bins for the dots}
}
\value{
A list of ggplot objects
}
\description{
Return individual plots for specified node for specimens within a GatingSet.
}
\examples{

library(flowCore)
library(flowWorkspace)
library(openCyto)
library(data.table)

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

SingleSpecimen <- Utility_IterativeGating(x=MyGatingSet[1],
 sample.name = "GUID", removestrings = removestrings,
 subset = "nonDebris", gate = "lymphocytes", xValue = "FSC-A",
  yValue = "SSC-A", bins = 270)

}
