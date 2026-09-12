#' Internal for Stacked Report
#' @importFrom dplyr select
#' @importFrom tidyselect all_of
#'
#' @noRd
ReportProcess <- function(x, columns){
  # x <- reports[[1]]
  data <- x %>% dplyr::select(all_of(columns))
}