#' Dashboard Internal, returns designated hex color.
#'
#' @param x The instrument designation
#' @param data The QC status data.frame
#'
#' @importFrom dplyr filter
#' @importFrom dplyr pull
#'
#' @return A hex code to fill with
#' @noRd
ColorCode <- function(x, data){

  Color <- data %>% dplyr::filter(Instrument %in% x) %>%
    pull(QCStatus)

  Hex <- "#FFFFFF"

  if (Color == "Red"){Hex <- "#C80815"}
  if (Color == "Orange"){Hex <- "#FF6E00"}
  if (Color == "Yellow"){Hex <- "#BA8E23"}
  if (Color == "Green"){Hex <- "#0B6623"}
  return(Hex)
}