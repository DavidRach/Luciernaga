#' Dashboard Internal, summarizes 3 months QC data for all instruments
#'
#' @param x A vector of instrument names
#' @param y A list of LevyJenningParse updated data objects
#' @param timewindow The number  desired months
#' 
#' @importFrom purrr map2
#' @importFrom dplyr bind_rows mutate filter pull group_by slice ungroup n
#' @importFrom tidyr pivot_wider
#'
#' @return Data ready for gt coloring
#' @noRd
QCHistory <- function(x, y, timewindow = 24) {
  InstrumentLength <- length(y)

  if (InstrumentLength > 1) {
    TheInstrumentLength <- 2
  } else {
    TheInstrumentLength <- 1
  }

  TheDataset <- map2(.x = x, .f = AcrossTime, .y = y) |> bind_rows()

  TheDates <- TheDataset |>
    group_by(DateTime) |>
    mutate(TheInstrumentCount = n()) |>
    filter(TheInstrumentCount >= TheInstrumentLength) |>
    pull(DateTime)

  Assembled <- TheDataset |> filter(DateTime %in% TheDates)

  Assembled <- Assembled |>
    group_by(DateTime, Instrument) |>
    slice(1) |>
    ungroup()

  Figure <- Assembled |>
    group_by(Instrument) |>
    pivot_wider(names_from = DateTime, values_from = QCStatus)

  return(Figure)
}