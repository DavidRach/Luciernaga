#' Wrapper for LuciernagaQC for large scale profiling. Matches unstained with single colors for respective
#' dates, extracts signature and plots, returns data.frame and list of patchwork objects.
#' 
#' @param x The list of Dates
#' @param UnstainedList The list of file.paths for the unstained files
#' @param FluorophoreList The list of file.paths for the fluorophore files
#' @param Multiple Default FALSE, if expecting multiple single color controls per unstained set TRUE
#' @param GateTemplatePath File.path to the openCyto gating template .csv
#' @param removestrings Default is ".fcs", removes from name
#' @param AFOverlap File.path to LuciernagaQC AFOverlap .csv to handle exceptions. 
#' @param controlType Either "beads" or "cells" (selects respective external AF protocol)
#' @param subsets The desired openCyto gating population to extract signature from
#' 
#' @importFrom purrr map
#' @importFrom dplyr bind_rows
#' @importFrom progressr with_progress
#' @importFrom progressr progressor
#' 
#' @export
#' 
#' @return A data.frame and a list of patchwork plots
#' 
#' @examples 
#' 
#' library(Luciernaga)
Luciernaga_SignatureExternalUnstained <- function(x, UnstainedList, FluorophoreList, Multiple=FALSE,
  GateTemplatePath, removestrings=".fcs", AFOverlap, controlType, subsets){
  with_progress({
    p <- progressor(along = x)  
    
    Iterated <- map(.x = x, .f = function(date) {

      result <- ControlWrapper(date, UnstainedList = UnstainedList, FluorophoreList = FluorophoreList,
                            Multiple = Multiple, GateTemplatePath = GateTemplatePath,
                            AFOverlap = AFOverlap, controlType=controlType, subsets=subsets)
      p()
      return(result)
    })
  })

  ThePlots <- lapply(Iterated, function(x) x$Plot)
  TheData <- lapply(Iterated, function(x) x$data)
  TheDataset <- bind_rows(TheData)
  TheFinalList <- list(data=TheDataset, plot=ThePlots)
  return(TheFinalList)
}

#' Internal for Luciernaga_SignatureExternalUnstained
#' 
#' @param x The list of Dates
#' @param UnstainedList The list of file.paths for the unstained files
#' @param FluorophoreList The list of file.paths for the fluorophore files
#' @param Multiple Default FALSE, if expecting multiple single color controls per unstained set TRUE
#' @param GateTemplatePath File.path to the openCyto gating template .csv
#' @param removestrings Default is ".fcs", removes from name
#' @param AFOverlap File.path to LuciernagaQC AFOverlap .csv to handle exceptions. 
#' @param controlType Either "beads" or "cells" (selects respective external AF protocol)
#' @param subsets The desired openCyto gating population to extract signature from
#' 
#' @importFrom stringr str_detect
#' @importFrom flowWorkspace load_cytoset_from_fcs
#' @importFrom flowWorkspace GatingSet
#' @importFrom data.table fread
#' @importFrom openCyto gatingTemplate
#' @importFrom openCyto gt_gating
#' @importFrom purrr map
#' @importFrom Biobase pData
#' @importFrom Biobase exprs
#' @importFrom BiocGenerics subset
#' @importFrom dplyr bind_rows
#' 
#' @return A list containing data and plots
#' 
#' @noRd
ControlWrapper <- function(x, UnstainedList, FluorophoreList, Multiple=FALSE,
  GateTemplatePath, removestrings=".fcs", AFOverlap, controlType, subsets){

  Value <- paste0(x, "(?!\\d)")
  Unstained <- UnstainedList[str_detect(UnstainedList, Value, negate = FALSE)]
  Fluorophore <- FluorophoreList[str_detect(FluorophoreList, Value, negate = FALSE)]

  if (Multiple == FALSE){
   if (length(Fluorophore) != length(Unstained)){
      stop("Mismatch unstained vs stained")}
  }

  files <- c(Unstained, Fluorophore)
  cs <- load_cytoset_from_fcs(files, truncate_max_range = FALSE, transformation = FALSE)
  gs <- GatingSet(cs)

  TheGates <- data.table::fread(GateTemplatePath)
  TheGating <- gatingTemplate(TheGates)
  gt_gating(TheGating, gs)

  Plot <- map(.x=gs, .f=Utility_GatingPlots, sample.name=c("GROUPNAME", "TUBENAME"),
   removestrings=removestrings, subset="root", gtFile=TheGates, DesiredGates=NULL,
  outpath = NULL, returnType="patchwork", therows=2, thecolumns=2)

  pd <- pData(gs)
  Names <- pData(gs)$name
  Unstained <- str_detect(Names, "Unstained")
  Unstained <- data.frame(Unstained, check.names=FALSE)
  newpd <- cbind(pd, Unstained)
  pData(gs) <- newpd

  TheUnstained <- subset(gs, Unstained == TRUE)
  TheFluorophore <- subset(gs, Unstained == FALSE)

  if (controlType == "beads"){

  UnstainedBeads <- Luciernaga_QC(x=TheUnstained[1], desiredAF = NULL,
                  subsets=subsets, removestrings=removestrings,
                  sample.name="GUID", unmixingcontroltype = "beads",
                  Unstained = TRUE, ratiopopcutoff = 0.01, Verbose = FALSE,
                  AFOverlap = AFOverlap, stats = "median",
                  ExportType = "data.frame", SignatureReturnNow = TRUE,
                  outpath = NULL)
  
  FluorophoreSignature <- map(.x=TheFluorophore, .f=Luciernaga_QC,
      subsets=subsets, removestrings=removestrings, sample.name="GUID",
      unmixingcontroltype = "beads", Unstained = FALSE, ratiopopcutoff = 0.01,
      Verbose = FALSE, AFOverlap = AFOverlap, stats = "median", 
      ExportType = "data", SignatureReturnNow = FALSE,
      outpath = NULL, Increments=0.1, SecondaryPeaks=2, experiment.name = "$DATE",
      condition = x, SCData="subtracted", NegativeType="default",
      BeadAF=UnstainedBeads, BeadMainAF="UV1-A") |> bind_rows()
  } else {
    UnstainedCells <- Luciernaga_QC(x=TheUnstained[1], desiredAF = NULL,
      subsets=subsets, removestrings=removestrings,
      sample.name="GUID", unmixingcontroltype = "cells",
      Unstained = TRUE, ratiopopcutoff = 0.01, Verbose = FALSE,
      AFOverlap = AFOverlap, stats = "median",
      ExportType = "data.frame", SignatureReturnNow = TRUE,
      outpath = NULL)

    FluorophoreSignature <- map(.x=TheFluorophore, .f=Luciernaga_QC,
      subsets=subsets, removestrings=removestrings, sample.name="GUID",
      unmixingcontroltype = "cells", Unstained = FALSE, ratiopopcutoff = 0.01,
      Verbose = FALSE, AFOverlap = AFOverlap, stats = "median", 
      ExportType = "data", SignatureReturnNow = FALSE,
      outpath = NULL, Increments=0.1, SecondaryPeaks=2, experiment.name = "$DATE",
      condition = x, SCData="subtracted", NegativeType="default",
      CellAF=UnstainedCells, CellMainAF="UV1-A") |> bind_rows()
  }
  
  TheReturns <- list(Plot=Plot, data=FluorophoreSignature)

  return(TheReturns)
}

