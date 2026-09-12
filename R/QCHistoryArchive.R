#' Dashboard Internal, summarizes x months QC data for all instruments
#'
#' @param x A vector of instrument names
#' @param y A list of LevyJenningParse updated data objects
#' @param timewindow The number  desired months
#' 
#' @importFrom purrr map
#' @importFrom dplyr bind_rows group_by mutate filter pull slice ungroup
#' @importFrom tidyr pivot_wider
#'
#' @return Data ready for gt coloring
#' @noRd
QCHistoryArchive <- function(x, historydata, timewindow=24){
  InstrumentLength <- length(x)

  if (InstrumentLength > 1){TheInstrumentLength <- 2
  } else {TheInstrumentLength <- 1}

  TheDataset <- map(.x=x, data=historydata,
     .f=InternalColorCodeStatus) |> bind_rows()

  TheDates <- TheDataset |> group_by(DateTime) |>
    mutate(TheInstrumentCount = n()) |>
    filter(TheInstrumentCount >= TheInstrumentLength) |> pull(DateTime)

  Assembled <- TheDataset |> filter(DateTime %in% TheDates)

  Assembled <- Assembled |> group_by(DateTime, Instrument) |> slice(1) |> ungroup()

  Figure <- Assembled |> group_by(Instrument) |>
    pivot_wider(names_from = DateTime, values_from = QCStatus)

  return(Figure)
}