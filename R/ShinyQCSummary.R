#' Dashboard Internal, processes to did parameter pass in past week
#'
#' @param x The data.frame output from LevyJennings or QCBeads Parse
#'
#' @importFrom dplyr mutate
#' @importFrom dplyr relocate
#' @importFrom dplyr pull
#' @importFrom purrr map
#' @importFrom dplyr bind_rows
#'
#' @return Data frame of passing status for respective parameters
#' @noRd
ShinyQCSummary <- function(x, Instrument){

  Intermediate <- x %>% mutate(Date=as.Date(DateTime)) %>%
    relocate(Date, .before=1)

  Dates <- Intermediate %>% pull(Date) %>% unique()

  Data <- map(.x=Dates, .f=ShinyQCSummaryParser, Intermediate=Intermediate) %>%
    bind_rows()

  Data <- Data %>% mutate(Instrument=Instrument) %>%
    relocate(Instrument, .after="Date")

  return(Data)
}