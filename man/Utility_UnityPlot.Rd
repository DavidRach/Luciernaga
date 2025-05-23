% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Utility_UnityPlots.R
\name{Utility_UnityPlot}
\alias{Utility_UnityPlot}
\title{Compare the plots for given parameters for all specimens in the gating set.}
\usage{
Utility_UnityPlot(
  x,
  y,
  GatingSet,
  marginsubset,
  gatesubset,
  sample.name,
  removestrings,
  clearance,
  bins,
  gatelines,
  reference,
  returntype,
  outpath,
  filename = NULL,
  cartesian = TRUE
)
}
\arguments{
\item{x}{The desired x-axis parameter}

\item{y}{The desired y-axis parameter}

\item{GatingSet}{The Gating Set that contains the specimens}

\item{marginsubset}{The desired gate for the margins}

\item{gatesubset}{The desired gate of cells you want to plot}

\item{sample.name}{The keyword that designates different specimens}

\item{removestrings}{Values to remove from both the plot titles and the pdf}

\item{clearance}{A value of clearance multiplied to the margin}

\item{bins}{How many bins to class the cells into}

\item{gatelines}{Whether to plot the reference lines}

\item{reference}{Reference for the gate cutoff lines}

\item{returntype}{Whether to return "pdf", "patchwork" or "plots"}

\item{outpath}{The desired location to send the assembled pdf to}

\item{filename}{Default NULL, provide name to set the filename.}

\item{cartesian}{Default TRUE, set to false to remove cartesian_coord centering}
}
\value{
The ggplots for all the specimens, as well as the optional .pdf
}
\description{
Compare the plots for given parameters for all specimens in the gating set.
}
\examples{

library(BiocGenerics)
library(flowCore)
library(flowWorkspace)
library(openCyto)
library(data.table)

File_Location <- system.file("extdata", package = "Luciernaga")
FCS_Files <- list.files(path = File_Location, pattern = ".fcs",
  full.names = TRUE)
Unmixed_FullStained <- FCS_Files[grep("Unmixed", FCS_Files)]
UnmixedFCSFiles <- Unmixed_FullStained[1:2]
UnmixedCytoSet <- load_cytoset_from_fcs(UnmixedFCSFiles[1:2],
  truncate_max_range = FALSE,transformation = FALSE)
UnmixedGatingSet <- GatingSet(UnmixedCytoSet)
Markers <- colnames(UnmixedCytoSet)
KeptMarkers <- Markers[-grep("Time|FS|SC|SS|Original|-W$|-H$|AF", Markers)]
MyBiexponentialTransform <- flowjo_biexp_trans(channelRange = 256,
  maxValue = 1000000,pos = 4.5, neg = 0, widthBasis = -1000)
TransformList <- transformerList(KeptMarkers, MyBiexponentialTransform)
UnmixedGatingSet <- flowWorkspace::transform(UnmixedGatingSet, TransformList)
FileLocation <- system.file("extdata", package = "Luciernaga")
UnmixedGates <- fread(file.path(path = FileLocation, pattern = 'GatesUnmixed.csv'))
UnmixedGating <- gatingTemplate(UnmixedGates)
gt_gating(UnmixedGating, UnmixedGatingSet)

removestrings <-  c("DTR_", ".fcs")
StorageLocation <- file.path("C:", "Users", "JohnDoe", "Desktop")

SingleUnityPlot <- Utility_UnityPlot(x="Spark Blue 550-A", y="BUV805-A",
GatingSet=UnmixedGatingSet, sample.name="GROUPNAME", bins=100, clearance=0.2,
removestrings=removestrings, marginsubset="lymphocytes", gatesubset="live",
gatelines=FALSE, reference=NULL, returntype="patchwork",outpath=StorageLocation)

}
