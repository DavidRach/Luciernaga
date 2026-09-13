#' Visualize cosine similarity of raw .fcs files to evaluate
#'  single color controls.
#'
#' @param path The location to the folder where the Luciernaga
#'  .fcs files are stored
#' @param reference A path to a .csv file or a dataframe containing
#'  Fluorophore and
#' Detector column information for the panel.
#' @param stats Whether to use the median or mean for fluorescent intensity.
#' @param LinePlots Return this kind of plot, default is set to TRUE
#' @param CosinePlots Return this kind of plot, default is set to TRUE
#' @param StackedBarPlots Return this kind of plot, default is set to TRUE
#' @param HeatmapPlots Return this kind of plot, default is set to TRUE
#' @param RetainedType Whether the data.frame contains "raw" or
#'  "normalized" values
#' @param experiment Provide directly experiment name (ex. "JAN2024")
#' @param condition Provide directly experiment name (ex. "JAN2024")
#' @param TheSummary Whether summarized (TRUE) or individual cells (FALSE).
#'
#' @importFrom dplyr select pull bind_rows mutate relocate
#' @importFrom purrr map
#' @importFrom tidyr separate
#' @importFrom utils read.csv
#'
#' @return A data.frame compatible with LuciernagaReport()
#'
#' @export
#'
#' @examples
#'
#' library(flowCore)
#' library(flowWorkspace)
#' library(openCyto)
#' library(data.table)
#' library(dplyr)
#' library(purrr)
#' library(stringr)
#'
#' StorageLocation <- file.path(tempdir(), "LuciernagaFCSToReportExample")
#' if (!dir.exists(StorageLocation)) {dir.create(StorageLocation)}
#'
#' File_Location <- system.file("extdata", package = "Luciernaga")
#' FCS_Files <- list.files(path = File_Location, pattern = ".fcs",
#'   full.names = TRUE)
#' CellSingleColorFiles <- FCS_Files[grep("Cells", FCS_Files)]
#' CellSingleColors <- CellSingleColorFiles[!str_detect("Unstained",
#'  CellSingleColorFiles)]
#' MyCytoSet <- load_cytoset_from_fcs(CellSingleColors[1:2],
#'   truncate_max_range = FALSE,transformation = FALSE)
#' MyGatingSet <- GatingSet(MyCytoSet)
#' MyGates <- fread(file.path(path = File_Location, pattern = 'Gates.csv'))
#' MyGatingTemplate <- gatingTemplate(MyGates)
#' gt_gating(MyGatingTemplate, MyGatingSet)
#' removestrings <-  c(".fcs", "(", ")", "Cells")
#'
#' FileLocation <- system.file("extdata", package = "Luciernaga")
#' pattern = "AutofluorescentOverlaps.csv"
#' AFOverlap <- list.files(path=FileLocation, pattern=pattern,
#'  full.names = TRUE)
#'
#' SingleColor_Data <- map(.x=MyGatingSet[1:2], .f=Luciernaga_QC,
#'  subsets="lymphocytes", removestrings=removestrings, sample.name="GUID",
#'  unmixingcontroltype = "cells", Unstained = FALSE, ratiopopcutoff = 0.001,
#'  Verbose = FALSE, AFOverlap = AFOverlap, stats = "median", ExportType = "fcs",
#'  Brightness=TRUE, SignatureReturnNow = FALSE,outpath = StorageLocation,
#'  Increments=0.1, SecondaryPeaks=2, experiment = "FirstExperiment",
#'  condition = "ILTPanel", Subtraction = "Internal", CellAF=TheCellAF,
#'   SCData="subtracted",NegativeType="default")
#'
#' TheLuciernagaOutputs_FCS <- list.files(StorageLocation, pattern="fcs",
#'  full.names = TRUE)
#' TheLuciernagaOutputs_CSV <- list.files(StorageLocation, pattern="csv",
#'  full.names = TRUE)
#' PanelPath <- file.path(File_Location, "Panel.csv")
#'
#' ReportOutput <- Luciernaga_FCSToReport(path=StorageLocation,
#'  reference=PanelPath,
#'  stats="median", RetainedType = "normalized", experiment="FirstExperiment",
#'  condition="ILTExperiment", TheSummary = TRUE)
#'
Luciernaga_FCSToReport <- function(path, reference, stats = "median",
                                    LinePlots = TRUE, CosinePlots = TRUE,
                                    StackedBarPlots = TRUE,
                                    HeatmapPlots = TRUE, RetainedType,
                                    experiment, condition,
                                    TheSummary = TRUE) {
  if (!is.data.frame(reference)) {
    CSV <- read.csv(reference, check.names = FALSE)
  } else {
    CSV <- reference
  }

  internalstrings <- c("-A")
  CSV$Fluorophore <- NameCleanUp(name = CSV$Fluorophore,
    removestrings = internalstrings)
  CSV$Detector <- NameCleanUp(name = CSV$Detector,
    removestrings = internalstrings)
  Variables <- CSV |> dplyr::select(Fluorophore) |> pull()
  fcsfiles <- list.files(path, pattern = ".fcs", full.names = TRUE)
  #x <- Variables[19]
  #inputfiles <- fcsfiles

  TheseFluorophores <- map(.x = Variables, .f = FluorophoreFilePresent,
    inputfiles = fcsfiles)
  TheseFluorophores <- Filter(Negate(is.null), TheseFluorophores)
  TheseFluorophores <- unlist(TheseFluorophores)
  #x <- TheseFluorophores[1]
  #data <- CSV
  #inputfiles = fcsfiles

  TheData <- map(.x = TheseFluorophores, .f = FCSImport, data = CSV,
    inputfiles = fcsfiles, RetainedType = RetainedType, stats = stats,
    TheSummary = TheSummary) |> bind_rows()

  TheData$Cluster <- gsub(" (Cells)", "", fixed = TRUE, TheData$Cluster)
  TheData$Cluster <- gsub(" (Beads)", "", fixed = TRUE, TheData$Cluster)

  if (TheSummary == TRUE) {
    TheExperiment <- as.character(experiment)
    TheCondition <- as.character(condition)

    TheData1 <- TheData |>
      separate(Cluster, into = c("Sample", "Cluster"), sep = "_")
    TheData1 <- TheData1 |>
      mutate(Experiment = TheExperiment)
    TheData1 <- TheData1 |>
      mutate(Condition = TheCondition)
    TheData1 <- TheData1 |>
      relocate(Sample, Experiment, Condition, .before = Cluster)

    TheData <- TheData1

    return(TheData)
  } else {
    TheData <- TheData |>
      separate(Cluster, into = c("Sample", "Cluster"), sep = "_")
    return(TheData)
  }
}