#' Converts the Luciernaga outputs into .pdf plots
#'
#' @param data The data.frame output from LuciernagaQC
#' @param RetainedType Whether the data.frame contains "raw" or
#'  "normalized" values
#' @param CellPopRatio What mininum ratio needed to retain cluster.
#' @param outfolder The location that you want to save the .pdf output to.
#' @param filename The name you want to save your .pdf file as.
#' @param LinePlots Passed to Utility_Patchwork for "pdf" or
#'  "patchwork" or "plots"
#' @param CosinePlots Return this kind of plot, default is set to TRUE
#' @param StackedBarPlots Return this kind of plot, default is set to TRUE
#' @param HeatmapPlots Return this kind of plot, default is set to TRUE
#' @param returntype Return "pdf", "patchwork" or "plots"
#' @param reference path or data.frame containing Fluorophore column
#'  for ordering
#' @param thecolumns The number of columns per page
#' @param therows The number of rows per page
#' @param width Desired page width
#' @param height Desired page height
#'
#' @importFrom dplyr group_by summarize left_join relocate
#'  filter mutate rename across select bind_rows bind_cols ungroup pull
#' @importFrom tidyselect all_of everything
#' @importFrom purrr map
#' @importFrom utils head tail read.csv
#'
#' @return A value to be determined later
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
#' File_Location <- system.file("extdata", package = "Luciernaga")
#' FCS_Files <- list.files(path = File_Location, pattern = ".fcs",
#'   full.names = TRUE)
#' CellSingleColorFiles <- FCS_Files[grep("Cells", FCS_Files)]
#' CellSingleColors <- CellSingleColorFiles[!str_detect("Unstained", CellSingleColorFiles)]
#' MyCytoSet <- load_cytoset_from_fcs(CellSingleColors[1:2],
#'   truncate_max_range = FALSE,transformation = FALSE)
#' MyGatingSet <- GatingSet(MyCytoSet)
#' MyGates <- fread(file.path(path = File_Location, pattern = 'Gates.csv'))
#' MyGatingTemplate <- gatingTemplate(MyGates)
#' gt_gating(MyGatingTemplate, MyGatingSet)
#' removestrings <-  c("DR_", "Cells", ".fcs", "-", " ")
#' StorageLocation <- file.path("C:", "Users", "JohnDoe", "Desktop")
#'
#' FileLocation <- system.file("extdata", package = "Luciernaga")
#' pattern = "AutofluorescentOverlaps.csv"
#' AFOverlap <- list.files(path=FileLocation, pattern=pattern, full.names = TRUE)
#'
#' SingleColor_Data <- map(.x=MyGatingSet[1:2], .f=Luciernaga_QC, subsets="lymphocytes",
#'  removestrings=removestrings, sample.name="GUID", unmixingcontroltype = "cells",
#'  Unstained = FALSE, ratiopopcutoff = 0.001, Verbose = FALSE, AFOverlap = AFOverlap,
#'  stats = "median", ExportType = "data", SignatureReturnNow = FALSE,
#'  outpath = TemporaryFolder, Increments=0.1, SecondaryPeaks=2,
#'  experiment = "FirstExperiment", condition = "ILTPanel", Subtraction = "Internal",
#'  CellAF=TheCellAF, SCData="subtracted", NegativeType="default") %>% bind_rows()
#'
#' pattern = "^Panel.csv"
#' CSV <- list.files(path=FileLocation, pattern=pattern, full.names=TRUE)
#' TheFluorophoreOrder <- read.csv(CSV, check.names = FALSE)
#'
#' ThePlots <- Luciernaga_Plots(data=SingleColor_Data, RetainedType="normalized",
#'  CellPopRatio=0.05, outfolder=NULL, filename="LuciernagaReport", returntype="plots",
#'  LinePlots=FALSE, CosinePlots=FALSE, StackedBarPlots = FALSE, HeatmapPlots = TRUE,
#'  reference = TheFluorophoreOrder)
#'
Luciernaga_Plots <- function(data, RetainedType, CellPopRatio, outfolder,
                              filename, LinePlots = TRUE, CosinePlots = TRUE,
                              StackedBarPlots = TRUE, HeatmapPlots = TRUE,
                              returntype = "patchwork", reference = NULL,
                              thecolumns = 2, therows = 2, width = 9,
                              height = 7) {
  if (!is.null(reference)) {
    if (!is.data.frame(reference)) {
      reference <- read.csv(reference, check.names = FALSE)
    }

    PreferredOrder <- reference |> pull(Fluorophore)
    PreferredOrder <- gsub("-A", "", PreferredOrder)
  } else {
    PreferredOrder <- NULL
  }

  #################################################
  # Filtered by CellPopRatio, and creating other  #
  #################################################

  TheCounts <- data |> group_by(Sample, Experiment, Condition) |>
    summarize(TotalCells = sum(Count, na.rm = TRUE), .groups = "drop")

  TheData <- data |>
    left_join(TheCounts, by = c("Sample", "Experiment", "Condition"))

  TheData <- TheData |> mutate(Ratio = round(Count / TotalCells, 3)) |>
    relocate(Ratio, .after = Count) |> select(-TotalCells)

  FilteredData <- TheData |> filter(Ratio > CellPopRatio)

  OtherData <- FilteredData |> group_by(Sample, Experiment, Condition) |>
    summarize(LostRatio = 1 - sum(Ratio, na.rm = TRUE), .groups = "drop")

  Other <- TheCounts |> left_join(OtherData, by = c("Sample", "Experiment",
    "Condition")) |> mutate(Count = round(TotalCells * LostRatio, 0)) |>
    select(-TotalCells) |> relocate(Count, .before = LostRatio) |>
    rename(Ratio = LostRatio) |> mutate(Cluster = "Other")

  OtherN <- nrow(Other)
  FirstDetectorColumn <- which(grepl("\\d", colnames(data)))[1]
  LastDetectorColumn <- tail(which(grepl("\\d", colnames(data))), 1)

  Replacement <- data |>
    select(all_of(FirstDetectorColumn:LastDetectorColumn)) |>
    head(OtherN) |> mutate(across(everything(), ~0))

  Replacements <- bind_cols(Other, Replacement) |> ungroup()
  Replaced <- bind_rows(FilteredData, Replacements)

  ##############
  # Lets Begin #
  ##############

  Items <- data.frame(table(data$Sample)) |> pull(Var1) |> as.character()

  if (!is.null(PreferredOrder)) {
    if (all(Items %in% PreferredOrder)) {
      Items <- PreferredOrder
    } else {
      message("names not matching, no reorderring according to panel order")
    }
  }

  #x <- Items[1]
  #data <- Replaced

  ThePlots <- map(.x = Items, .f = Luciernaga:::InternalReport, data = Replaced,
    FirstDetectorColumn = FirstDetectorColumn,
    LastDetectorColumn = LastDetectorColumn,
    RetainedType = RetainedType, CellPopRatio = CellPopRatio,
    LinePlots = LinePlots, CosinePlots = CosinePlots,
    StackedBarPlots = StackedBarPlots, HeatmapPlots = HeatmapPlots)

  if (returntype == "pdf") {
    Utility_Patchwork(x = ThePlots, filename = filename, outfolder = outfolder,
      thecolumns = thecolumns, therows = therows, width = width,
      height = height, returntype = "pdf", NotListofList = FALSE)
  }

  if (returntype == "patchwork") {
    Hey <- Utility_Patchwork(x = ThePlots, filename = filename,
      outfolder = outfolder, thecolumns = thecolumns, therows = therows,
      width = width, height = height, returntype = "patchwork",
      NotListofList = FALSE)
    return(Hey)
  }

  if (returntype == "plots") {
    return(ThePlots)
  }
}