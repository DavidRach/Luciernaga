#' Small Internal Function
#' 
#' @param x Something
#' 
#' @importFrom dplyr case_when
#' 
#' @return Something
#' 
#' @noRd                       
InstrumentText <- function(x) {
  dplyr::case_when(
    x == "Green" ~ "Pass",
    x == "Yellow" ~ "Caution",
    x == "Orange" ~ "Caution",
    x == "Red" ~ "Fail",
    TRUE ~ NA_character_)
}