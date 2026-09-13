#' Given either Instruments or Fluors, compares cosine values
#'  across instruments
#' for all or selected fluorophores against the provided MainFluorophore
#'
#' @param MainFluorophore The desired fluorophore to compare against
#' @param Instruments A list of instrument detectors for the 
#' respective instruments
#' @param TheseFluorophores A vector of fluorophore names to compare,
#'  when Instruments = NULL
#' @param returnType Whether to return "data" or alternatively a "plot"
#'
#' @importFrom purrr map
#' @importFrom dplyr bind_rows filter group_by arrange desc
#'  select slice ungroup pull bind_cols
#' @importFrom tidyr pivot_wider
#'
#' @return A data.frame or in the future a plot object
#'
#' @examples
#'
#' library(dplyr)
#' library(tidyr)
#' library(purrr)
#'
#' Data <- Luciernaga:::InstrumentComparison(
#' TheseFluorophores=c("BV421", "BV510", "BV605", "BV650"), 
#' MainFluorophore="FITC", returnType="data")
#'
#' @noRd
InstrumentComparison <- function(MainFluorophore,
                                  Instruments = NULL,
                                  TheseFluorophores = NULL,
                                  returnType) {

  if (is.null(Instruments) && is.null(TheseFluorophores)) {
    stop("Provide either an Instruments or TheseFluorophores argument")
  }
  if (!is.null(Instruments) && !is.null(TheseFluorophores)) {
    message("Typically select data using either Instrument or
     TheseFluorophores, set other to NULL. Will prioritize
     selecting TheseFluorophores and only show Instruments
     for which there are references")
  }

  if (is.null(TheseFluorophores)) {
    Values <- map(.x = Instruments, .f = Luciernaga:::InstrumentNameCheck)
    Hmm <- unlist(Values)
    TheFluorophores <- map(.x = Hmm, .f = Luciernaga:::InstrumentReturn)
    Shared <- Reduce(intersect, TheFluorophores)
  } else {
    Hmm <- GetInstruments()
    TheReferences <- map(.x = Hmm, .f = Luciernaga:::InstrumentReturn)
    CheckThese <- c(MainFluorophore, TheseFluorophores)
    LogicalReturn <- sapply(TheReferences,
                             function(ref) all(CheckThese %in% ref))
    Hmm <- Hmm[LogicalReturn]
    Present <- TheReferences[LogicalReturn]
    Shared <- CheckThese
  }

  Dataset <- map(.x = Hmm, .f = Luciernaga:::InstrumentData,
                 fluorophores = Shared) |>
    bind_rows()
  # Dataset |> pull(Instrument) |> unique()

  if (!MainFluorophore %in% Shared) {
    stop(MainFluorophore, " not found across selected instruments")
  }

  if (!is.null(Instruments)) {
    OrderReference <- Hmm[1]
  } else {
    OrderReference <- Hmm[1]
  }

  ReferenceData <- Luciernaga:::InstrumentReferences(
    NumberDetectors = OrderReference)
  ReferenceLevels <- ReferenceData |>
    filter(Fluorophore %in% MainFluorophore) |>
    pull(Detector)
  ReferenceData$Detector <- factor(ReferenceData$Detector,
                                    levels = ReferenceLevels)

  TheReferenceList <- ReferenceData |> filter(Fluorophore %in% Shared)
  Shared_Rearranged <- TheReferenceList |>
    group_by(Fluorophore) |>
    arrange(desc(AdjustedY)) |>
    slice(1) |>
    select(Fluorophore, Detector) |>
    ungroup() |>
    arrange(Detector) |>
    pull(Fluorophore)

  LongReferences <- Dataset |>
    pivot_wider(names_from = "Detector", values_from = "AdjustedY")

  instrument_match <- LongReferences |> pull(Instrument) |> unique()

  Forward <- map(.x = instrument_match, .f = MatchRearrange,
                 thematch = Shared_Rearranged, data = LongReferences) |>
    bind_rows()
  # Forward |> pull(Instrument) |> unique()

  Comparison <- map(.x = instrument_match, .f = CosineReturn,
                     MainFluorophore = MainFluorophore, data = Forward) |>
    bind_cols()

  if (returnType == "data") {
    return(Comparison)
  } else if (returnType == "plot") {
    message("Add gt or ggplot2 option here")
    return(Comparison)
  } else {
    return(Comparison)
  }
}