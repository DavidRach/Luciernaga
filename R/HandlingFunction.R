#' Processes experiment folder for Luciernaga Signatures, intended
#' for David's panel, not yet generalized
#' 
#' @param x File.path to a target folder
#' 
#' @importFrom stringr str_detect
#' @importFrom data.table fread
#' @importFrom openCyto gatingTemplate
#' @importFrom flowWorkspace load_cytoset_from_fcs
#' @importFrom flowWorkspace GatingSet
#' @importFrom openCyto gt_gating
#' 
#' @return Luciernaga_QC processed signatures to a new Luciernaga subfolder
#' 
#' @noRd
#' 
#' @examples 
#' A <- 2+2
HandlingFunction <- function(x){

  # Isolate unmixing control types
  files <- list.files(x, pattern=".fcs", full.names=TRUE)
  Unstained <- files[grep("Unstained", files)]
  BeadsUnstained <- Unstained[grep("Beads", Unstained )]
  CellsUnstained <- Unstained[!stringr::str_detect(
      Unstained, "Beads|Dead")]

  Stained <- files[!stringr::str_detect(files, "Unstained")]
  BeadsStained <- Stained[grep("Beads", Stained)]
  CellsStained <- Stained[!stringr::str_detect(Stained, "Beads")]

  # References
  FileLocation <- system.file("extdata", package = "Luciernaga")
  CellGates <- data.table::fread(file.path(
      path = FileLocation, pattern = 'Gates.csv'))
  CellTemplate <- gatingTemplate(CellGates)
  BeadGates <- data.table::fread(file.path(
      path = FileLocation, pattern = 'GatesBeads.csv'))
  BeadTemplate <- gatingTemplate(BeadGates)
  pattern = "AutofluorescentOverlaps.csv"
  AFOverlap <- list.files(path=FileLocation,
   pattern=pattern, full.names = TRUE)

  # Gating
  UnstainedCell_CS <- load_cytoset_from_fcs(CellsUnstained,
   truncate_max_range = FALSE, transform = FALSE)
  UnstainedCell_GS <- GatingSet(UnstainedCell_CS)
  gt_gating(CellTemplate, UnstainedCell_GS)
  UnstainedCellPlot <- Utility_GatingPlots(x=UnstainedCell_GS[[1]],
   sample.name = "GUID", removestrings = ".fcs",
   gtFile = CellGates, DesiredGates = NULL,
   outpath = NULL, returnType="patchwork")

  StainedCell_CS <- load_cytoset_from_fcs(CellsStained,
   truncate_max_range = FALSE, transform = FALSE)
  StainedCell_GS <- GatingSet(StainedCell_CS)
  gt_gating(CellTemplate, StainedCell_GS)
  StainedCellPlot <- Utility_GatingPlots(x=StainedCell_GS[[1]],
   sample.name = "GUID", removestrings = ".fcs",
   gtFile = CellGates, DesiredGates = NULL,
   outpath = NULL, returnType="patchwork")

  UnstainedBeads_CS <- load_cytoset_from_fcs(BeadsUnstained,
   truncate_max_range = FALSE, transform = FALSE)
  UnstainedBeads_GS <- GatingSet(UnstainedBeads_CS)
  gt_gating(BeadTemplate, UnstainedBeads_GS)
  UnstainedBeadPlot <- Utility_GatingPlots(x=UnstainedBeads_GS[[1]],
   sample.name = "GUID", removestrings = ".fcs",
   gtFile = BeadGates, DesiredGates = NULL,
   outpath = NULL, returnType="patchwork")

  StainedBeads_CS <- load_cytoset_from_fcs(BeadsStained,
   truncate_max_range = FALSE, transform = FALSE)
  StainedBeads_GS <- GatingSet(StainedBeads_CS)
  gt_gating(BeadTemplate, StainedBeads_GS)
  StainedBeadPlot <- Utility_GatingPlots(x=StainedBeads_GS[[1]],
   sample.name = "GUID", removestrings = ".fcs",
   gtFile = BeadGates, DesiredGates = NULL,
   outpath = NULL, returnType="patchwork")

  These <- c(UnstainedCellPlot, StainedCellPlot, UnstainedBeadPlot, StainedBeadPlot)

  OnLocation <- list.dirs(x, full.names=TRUE, recursive=TRUE)
  Presence <- list.files(OnLocation, pattern="Luciernaga", full.names=TRUE)
  if (length(Presence) == 0){dir.create(file.path(x, "Luciernaga"))}
  Presence <- list.files(OnLocation, pattern="Luciernaga", full.names=TRUE)

  Utility_Patchwork(These, filename="GatingReport", outfolder=Presence,
  thecolumns=1, therows=2)

  # Single-color Cells

  UnstainedCellAF <- Luciernaga_QC(x=UnstainedCell_GS[[1]],
      subsets="nonDebris", removestrings=".fcs", sample.name="GUID",unmixingcontroltype = "cells", Unstained = FALSE, 
      ratiopopcutoff = 0.001, Verbose = FALSE, AFOverlap = AFOverlap,
      stats = "median", ExportType = "fcs", SignatureReturnNow = TRUE,
      outpath = Presence, Increments=0.1, minimalfcscutoff=0.001,
      experiment="SingleColor",condition="Cells")

  Unstained <- purrr::map(.f=Luciernaga_QC, .x=UnstainedCell_GS,
      subsets="nonDebris", removestrings=".fcs", sample.name="GUID",unmixingcontroltype = "cells", Unstained = FALSE, 
      ratiopopcutoff = 0.001, Verbose = FALSE, AFOverlap = AFOverlap,
      stats = "median", ExportType = "fcs", SignatureReturnNow = FALSE,
      outpath = Presence, Increments=0.1, minimalfcscutoff=0.001,
      experiment="SingleColor",condition="Cells")

  Colors <- purrr::map(.f=Luciernaga_QC, .x=StainedCell_GS,
      subsets="nonDebris", removestrings=".fcs", sample.name="GUID",unmixingcontroltype = "cells", Unstained = FALSE, 
      ratiopopcutoff = 0.001, Verbose = FALSE, AFOverlap = AFOverlap,
      stats = "median", ExportType = "fcs", SignatureReturnNow = FALSE,
      outpath = Presence, Increments=0.1, minimalfcscutoff=0.001,
      experiment="SingleColor",condition="Cells", CellAF=UnstainedCellAF)
}