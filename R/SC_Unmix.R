#' Scrambled Egg Unmixing of Single-Colors
#'
#' @param x A Gating Set object
#' @param sample.name The keyword containing the fcs file name
#' @param removestrings A list of values to remove from name
#' @param subset A gating hierarchy level to sort cells at, expression values retrieved
#' from these
#' @param multiplier A number to scale the OLS coefficients by
#' @param outpath The return folder for the .fcs files
#' @param returntype Whether to return "fcs" or "flowframe"
#' @param Verbose For troubleshooting name after removestrings
#' @param addon Additional addon to append to the new .fcs file name
#' @param ratiopopcutoff Desired cutoff for detector detection
#' @param NumberFluors Desired number additional fluors in matrix
#'
#' @importFrom flowCore keyword exprs write.FCS
#' @importFrom flowWorkspace gs_pop_get_data
#' @importFrom dplyr mutate select pull arrange filter
#' @importFrom tidyselect all_of where
#' @importFrom utils read.csv
#' @importFrom stats lsfit quantile
#' @importFrom purrr map
#' @importFrom rlang .data
#'
#' @return A new .fcs file with the new columns appended
#' @noRd
SC_Unmix <- function(x,
                      sample.name,
                      removestrings,
                      subset,
                      multiplier,
                      outpath,
                      returntype,
                      Verbose,
                      addon,
                      ratiopopcutoff,
                      NumberFluors = 1) {

  # Retrieving Single Color Sample data for Unmixing
  if (length(sample.name) == 2) {
    first <- sample.name[[1]]
    second <- sample.name[[2]]
    first <- keyword(x, first)
    second <- keyword(x, second)
    name <- paste(first, second, sep = "_")
  } else {
    name <- keyword(x, sample.name)
  }

  name <- NameCleanUp(name, removestrings = removestrings)
  if (Verbose == TRUE) {
    message("After removestrings, name is ", name)
  }

  cs <- gs_pop_get_data(x, subset)
  Data <- exprs(cs[[1]])
  Data <- data.frame(Data, check.names = FALSE)

  OriginalColumnsVector <- colnames(Data)
  OriginalColumns <- colnames(Data)
  OriginalColumns <- data.frame(OriginalColumns, check.names = FALSE)
  OriginalColumnsIndex <- OriginalColumns |>
    mutate(IndexLocation = 1:nrow(.)) # TODO: `.` refers to lhs (magrittr-only)

  Backups <- Data |>
    mutate(Backups = 1:nrow(Data)) |>
    select(Backups)

  StashedDF <- Data[, grep("Time|FS|SC|SS|Original|W$|H$", names(Data))]
  StashedDF <- cbind(Backups, StashedDF)

  TheSampleData <- Data[, -grep("Time|FS|SC|SS|Original|W$|H$", names(Data))]
  BackupNames <- colnames(TheSampleData)
  startingcells <- nrow(TheSampleData)
  TheColsN <- ncol(TheSampleData)

  # Retrieving Single Color Average Brightness
  PeakReturn <- Luciernaga:::PeakDetectors(NewData = TheSampleData)

  WorkAround <- PeakReturn[[1]]
  PeakDetectorCounts <- PeakReturn[[2]]
  CellCutoff <- startingcells * ratiopopcutoff
  Detectors <- PeakDetectorCounts |> filter(Counts > CellCutoff)
  Retained <- Luciernaga:::RetainTheDetectors(AFOverlap = AFOverlap,
                                                Detectors = Detectors,
                                                name = name)

  if (!length(Retained) > 1) {
    RetainedA <- paste0(Retained, "-A")
    MaxValues <- WorkAround |>
      filter(.data[[Retained]] == 1) |>
      select(all_of(RetainedA)) |>
      pull(.) # TODO: `.` refers to lhs (magrittr-only)
  } else {
    MaxValues <- map(.x = Retained, .f = RetainedParse, data = WorkAround) |>
      unlist()
  }
  # hist(MaxValues)
  TargetValue <- quantile(MaxValues, 0.95)
  TargetValue <- TargetValue[[1]]

  # Returning SC Data from References

  controlData <- Luciernaga:::ReferenceScramble(name = name,
                                                  NumberDetectors = TheColsN,
                                                  NumberFluors = NumberFluors)
  NewNames <- c("Fluorophore", colnames(TheSampleData))
  colnames(controlData) <- NewNames

  TheNames <- controlData |> select(!where(is.numeric))
  Numerics <- controlData |> select(where(is.numeric))
  Numerics <- Numerics * TargetValue
  MFIedControlData <- cbind(TheNames, Numerics)

  NewNames <- MFIedControlData |> select(Fluorophore) |> pull()
  TheControlData <- MFIedControlData |> select(!Fluorophore)

  LeastSquares <- lsfit(x = t(TheControlData), y = t(TheSampleData),
                         intercept = FALSE)
  UnmixedData <- t(LeastSquares$coefficients)
  UnmixedData2 <- UnmixedData * multiplier

  colnames(UnmixedData2) <- NewNames
  TheData <- cbind(StashedDF, UnmixedData2)
  TheData <- TheData |> select(-Backups)
  rownames(TheData) <- NULL

  Ligands <- NewNames

  new_fcs <- Luciernaga:::InternalUnmix(cs = cs, StashedDF = StashedDF,
                                          TheData = TheData, Ligands = Ligands)
  # View(new_fcs@description)

  if (!is.null(addon)) {
    name <- paste0(name, addon)
  }

  CurrentTime <- format(Sys.time(), "%H:%M:%S")
  CurrentTime <- gsub(":", "", CurrentTime)

  AssembledName <- paste0(name, "_", CurrentTime, ".fcs")

  new_fcs@description$GUID <- AssembledName
  new_fcs@description$`$FIL` <- AssembledName

  if (is.null(outpath)) {
    outpath <- getwd()
  }

  fileSpot <- file.path(outpath, AssembledName)

  if (returntype == "fcs") {
    write.FCS(new_fcs, filename = fileSpot, delimiter = "#")
  } else {
    return(new_fcs)
  }
}