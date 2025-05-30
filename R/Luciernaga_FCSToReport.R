#' Visualize cosine similarity of raw .fcs files to evaluate single color
#' controls.
#'
#' @param path The location to the folder where the Luciernaga .fcs files are stored
#' @param reference A path to a .csv file or a dataframe containing Fluorophore and
#' Detector column information for the panel.
#' @param stats Whether to use the median or mean for fluorescent intensity.
#' @param LinePlots Return this kind of plot, default is set to TRUE
#' @param CosinePlots Return this kind of plot, default is set to TRUE
#' @param StackedBarPlots Return this kind of plot, default is set to TRUE
#' @param HeatmapPlots Return this kind of plot, default is set to TRUE
#' @param RetainedType Whether the data.frame contains "raw" or "normalized" values
#' @param TheSummary Whether summarized (TRUE) or individual cells (FALSE).
#' @param experiment Provide directly experiment name (ex. "JAN2024")
#' @param condition Provide directly experiment name (ex. "JAN2024")
#'
#' @importFrom dplyr select
#' @importFrom dplyr pull
#' @importFrom purrr map
#' @importFrom dplyr bind_rows
#' @importFrom dplyr mutate
#' @importFrom dplyr relocate
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
#' CellSingleColors <- CellSingleColorFiles[!str_detect("Unstained", CellSingleColorFiles)]
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
#' AFOverlap <- list.files(path=FileLocation, pattern=pattern, full.names = TRUE)
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
#' TheLuciernagaOutputs_FCS <- list.files(StorageLocation, pattern="fcs", full.names = TRUE)
#' TheLuciernagaOutputs_CSV <- list.files(StorageLocation, pattern="csv", full.names = TRUE)
#' PanelPath <- file.path(File_Location, "Panel.csv")
#'
#' ReportOutput <- Luciernaga_FCSToReport(path=StorageLocation, reference=PanelPath,
#'  stats="median", RetainedType = "normalized", experiment="FirstExperiment",
#'  condition="ILTExperiment", TheSummary = TRUE)
#'
Luciernaga_FCSToReport <- function(path, reference, stats = "median",
                                    LinePlots = TRUE, CosinePlots = TRUE,
                                    StackedBarPlots = TRUE, HeatmapPlots = TRUE,
                                    RetainedType, experiment, condition,
                                    TheSummary = TRUE){

  if (!is.data.frame(reference)){CSV <- read.csv(reference, check.names = FALSE)
  } else {CSV <- reference}

  internalstrings <- c("-A")
  CSV$Fluorophore <- NameCleanUp(name=CSV$Fluorophore, removestrings=internalstrings)
  CSV$Detector <- NameCleanUp(name=CSV$Detector, removestrings=internalstrings)
  Variables <- CSV %>% dplyr::select(Fluorophore) %>% pull(.)
  fcsfiles <- list.files(path, pattern=".fcs", full.names = TRUE)
  #x <- Variables[19]
  #inputfiles <- fcsfiles

  TheseFluorophores <- map(.x=Variables, .f=FluorophoreFilePresent,
                           inputfiles = fcsfiles)
  TheseFluorophores <- Filter(Negate(is.null), TheseFluorophores)
  TheseFluorophores <- unlist(TheseFluorophores)
  #x <- TheseFluorophores[1]
  #data <- CSV
  #inputfiles = fcsfiles

  TheData <- map(.x = TheseFluorophores, .f = FCSImport, data = CSV,
                   inputfiles = fcsfiles, RetainedType=RetainedType, stats=stats,
                 TheSummary=TheSummary) %>% bind_rows()

  TheData$Cluster <- gsub(" (Cells)", "", fixed =TRUE, TheData$Cluster)
  TheData$Cluster <- gsub(" (Beads)", "", fixed =TRUE, TheData$Cluster)

  if (TheSummary == TRUE){
  TheExperiment <- as.character(experiment)
  TheCondition <- as.character(condition)

  TheData1 <- TheData %>% separate(Cluster, into = c("Sample", "Cluster"), sep = "_")
  TheData1 <- TheData1 %>% mutate(Experiment=TheExperiment)
  TheData1 <- TheData1 %>% mutate(Condition=TheCondition)
  TheData1 <- TheData1 %>% relocate(Sample, Experiment, Condition, .before=Cluster)

  TheData <- TheData1

  return(TheData)
  } else {TheData <- TheData %>%
    separate(Cluster, into = c("Sample", "Cluster"), sep = "_")
    return(TheData)}
}


#' Internal for LuciernagaReportFromFCS
#'
#' @param x Passed Fluorophore Name
#' @param inputfiles List of .fcs files from path
#'
#' @importFrom stringr str_detect
#'
#' @return An internal value
#'
#' @noRd
FluorophoreFilePresent <- function(x, inputfiles){
  fcs_files <- inputfiles[str_detect(basename(inputfiles), x) &
                            str_detect(basename(inputfiles), ".fcs$")]
  if (x %in% c("PE", "APC")){
      x <- paste0(x, "-")
      fcs_files <- fcs_files[!str_detect(basename(fcs_files), x)]
      x <- gsub("-", "", x)
    } #ExceptionHandling

  if (length(fcs_files) > 0){return(x)}
}


#' Internal for LuciernagaReportFromFCS
#'
#' @param x A passed single cytoset object
#' @param Fluorophore The detector
#'
#' @importFrom flowCore keyword
#' @importFrom flowCore exprs
#' @importFrom dplyr mutate
#' @importFrom dplyr relocate
#'
#' @return An internal value
#'
#' @noRd
FCSImportFile <- function(x, Fluorophore, sample.name = "FILENAME"){

  filename <- keyword(x, sample.name)
  filename <- sub(".*\\\\", "", filename)
  filename <- sub(paste0(".*", Fluorophore), Fluorophore, filename)
  filename <- gsub(".fcs$", "", filename)
  rownames(filename) <- NULL

  #df <- exprs(x[[1]])
  df <- exprs(x)
  TheDF <- data.frame(df, check.names = FALSE)
  TheDF <- TheDF[,-grep("Time|FS|SC|SS|Original|W$|H$", names(TheDF))]
  colnames(TheDF) <- gsub("-A$", "", colnames(TheDF))
  DFNames <- TheDF %>% mutate(Cluster = filename) %>% relocate(Cluster,
                                                               .before = 1)
  return(DFNames)
}

#' Internal for LuciernagaReportFromFCS
#'
#' @param x Passed Fluorophore Name
#' @param data Passed data.frame of Fluorophore with Detectors
#' @param inputfiles List of .fcs files from path
#' @param RetainedType Whether the data.frame contains "raw" or "normalized" values
#' @param TheSummary Whether to return individual cells or summarized by stats.
#' RetainedType Whether the data.frame contains "raw" or "normalized" values
#'
#' @importFrom dplyr filter
#' @importFrom dplyr pull
#' @importFrom stringr str_detect
#' @importFrom flowWorkspace load_cytoset_from_fcs
#' @importFrom purrr map
#' @importFrom dplyr mutate
#' @importFrom dplyr across
#' @importFrom tidyselect where
#' @importFrom dplyr select
#' @importFrom dplyr group_by
#' @importFrom dplyr ungroup
#' @importFrom dplyr arrange
#' @importFrom dplyr left_join
#'
#' @return An internal value
#'
#' @noRd
FCSImport <- function(x, data, inputfiles, RetainedType, TheSummary, stats){

  # For each Fluorophore
  TheDetector <- data %>% dplyr::filter(Fluorophore %in% x) %>% pull(Detector)

  fcs_files <- inputfiles[str_detect(basename(inputfiles), x) &
                            str_detect(basename(inputfiles), ".fcs$")]

  if (x %in% c("PE", "APC")){x <- paste0(x, "-") #ExceptionHandling
  fcs_files <- fcs_files[!str_detect(basename(fcs_files), x)]
  x <- gsub("-", "", x)
  } #ExceptionHandling

  if (!length(fcs_files) == 0){

    cs <- load_cytoset_from_fcs(fcs_files, truncate_max_range = FALSE,
                                transformation = FALSE)
    thex <- x

    # Retrieve exprs data from each .fcs file and create cluster column
    # x <- cs[1]
    TheData <- map(.x = cs, .f = FCSImportFile, Fluorophore = thex) %>%
      bind_rows()

    # Removing Any Artificial Negatives Inserted By Luciernaga
    TheData <- TheData %>%
      mutate(Summed = rowSums(across(where(is.numeric)), na.rm = TRUE))

    TheData <- TheData %>% group_by(Summed) %>% dplyr::filter(n() <= 5) %>%
      ungroup() %>% dplyr::filter(!Summed == 0) %>% dplyr::select(-Summed)

    # Return Summarized Data
    if (RetainedType == "normalized"){
      Cluster <- TheData %>% dplyr::select(Cluster)
      DetectorData <- TheData %>% dplyr::select(-Cluster)
      DetectorData[DetectorData < 0] <- 0
      A <- do.call(pmax, DetectorData)
      Normalized <- DetectorData/A
      Normalized <- round(Normalized, 3)
      TheData <- cbind(Cluster, Normalized)
    }

    if (TheSummary == TRUE){
      TheTable <- data.frame(table(TheData$Cluster), check.names = FALSE)
      TheTable <- TheTable %>% dplyr::arrange(desc(Freq))
      colnames(TheTable)[1] <- "Cluster"
      colnames(TheTable)[2] <- "Count"

      TheClusters <- TheTable$Cluster

      SmallHelper <- function(x, data, stats){
      Cluster <- x
      Internal <- data %>% dplyr::filter(Cluster %in% x) %>% dplyr::select(where(is.numeric))
      Summarized <- AveragedSignature(x=Internal, stats=stats)
      Summarized <- cbind(Cluster, Summarized)
      }

      Summarized <- map(.x=TheClusters, .f=SmallHelper, data = TheData,
                        stats=stats) %>% bind_rows()

      ReturnFrame <- left_join(TheTable, Summarized, by="Cluster")
      return(ReturnFrame)
    } else {return(TheData)}
  }

}
