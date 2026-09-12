#' Internal for QC_Retrieval
#'
#' @param x Passed argument 1
#' @param y Passed argument 2
#' @param TheData The datset
#'
#' @importFrom dplyr select
#' @importFrom tidyselect all_of
#' @importFrom dplyr slice
#'
#' @return An internal value
#'
#' @keywords internal
RetrievalMerge <- function(x, y, TheData){
  Individual <- TheData %>% select(all_of(c(x, y))) %>% slice(1)
  Cell <- Individual %>% pivot_wider(names_from = 1, values_from = 2)
  return(Cell)
}