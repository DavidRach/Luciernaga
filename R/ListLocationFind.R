#' Internal for CytoSet check, actual filtering of subset list
#' 
#' @param x The iterated identity
#' @param data Intermediate data.frame containing RowNumber column
#' @param TheList The original list of fcs files to be sorted from
#' 
#' @importFrom dplyr filter
#' @importFrom dplyr pull
#' 
#' @return A compatible subset list
#' 
#' @noRd
ListLocationFind <- function(x, data, TheList){
  TheRows <- data |> filter(ID %in% x) |> pull(Iteration)
  TheSubset <- TheList[TheRows]
  return(TheSubset)
}