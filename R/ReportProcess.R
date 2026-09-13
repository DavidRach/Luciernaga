#' Internal for Stacked Report
#'
#' @param x TBD
#' @param columns TBD
#'
#' @importFrom dplyr select
#' @importFrom tidyselect all_of
#'
#' @return A data.frame containing only the specified columns
#' @noRd
ReportProcess <- function(x, columns) {
  # x <- reports[[1]]
  data <- x |> select(all_of(columns))
}