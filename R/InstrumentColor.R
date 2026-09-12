#' Small Internal Function
#' 
#' @param x Something
#' 
#' @importFrom dplyr case_when
#' 
#' @return Something
#' 
#' @noRd
InstrumentColor <- function(x) {
  dplyr::case_when(
    x == "Green" ~ "success",
    x == "Yellow" ~ "caution",
    x == "Orange" ~ "warning",
    x == "Red" ~ "danger",
    TRUE ~ NA_character_)
}