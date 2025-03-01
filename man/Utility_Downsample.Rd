% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Utility_Downsample.R
\name{Utility_Downsample}
\alias{Utility_Downsample}
\title{Downsample a gs based on subset}
\usage{
Utility_Downsample(
  x,
  sample.name,
  removestrings,
  subsets,
  subsample = NULL,
  inverse.transform,
  internal = FALSE,
  export = FALSE,
  outpath = NULL,
  metadataCols = NULL
)
}
\arguments{
\item{x}{A gating set object}

\item{sample.name}{Keyword specifying sample name}

\item{removestrings}{Value to be removed from sample name}

\item{subsets}{The gating hierarchy subset you want to include}

\item{subsample}{Total number of events to sub-sample from each specimen.
If value between 0 and 1, grabs that equivalent proportion of total cells.}

\item{inverse.transform}{Whether to reverse the GatingSet Transform on the data,
default is set to FALSE.}

\item{internal}{Whether to return as a data.frame (vs. a flow frame or .fcs file)}

\item{export}{Default is set to FALSE, when TRUE (and internal = FALSE) returns a
.fcs file to outpath.}

\item{outpath}{When export is true, the file.path to where you want the .fcs file
stored.}

\item{metadataCols}{column names from pData to append as metadata for the .fcs, default NULL}
}
\value{
Either a data.frame, a flow.frame or an .fcs file depending on your
selected options
}
\description{
Downsample a gs based on subset
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

removestrings <- c("DTR_2023_ILT_15_Tetramers-","-Ctrl_Tetramer_Unmixed", ".fcs")
StorageLocation <- file.path("C:", "Users", "JohnDoe", "Desktop")

CountData <- gs_pop_get_count_fast(UnmixedGatingSet)
Counts_Specimen <- CountData \%>\%
 filter(Population \%in\% "/singletsFSC/singletsSSC/singletsSSCB/nonDebris/lymphocytes/live") \%>\%
 select(name, Count)

SingleSample <- Utility_Downsample(UnmixedGatingSet[1],
 sample.name = "GROUPNAME", removestrings=removestrings,
 subsets = "live", subsample = 2500, internal = FALSE,
 export = FALSE, inverse.transform=TRUE)

}
