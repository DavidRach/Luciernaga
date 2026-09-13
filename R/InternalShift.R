#' Internal for Column Append, ensures nothing is zero valued
#'
#' @param x The column name to select and shift
#' @param columnframe The data.frame containing the column
#'
#' @importFrom dplyr select
#' @importFrom tidyselect all_of
#'
#' @return An internal value
#'
#' @noRd
InternalShift <- function(x, columnframe) {

  TheColumn <- columnframe |> select(all_of(x))
  ShiftedColumn <- TheColumn + abs(min(TheColumn))+1
  return(ShiftedColumn)
}