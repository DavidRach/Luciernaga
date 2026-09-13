#' Internal for InstrumentComparison, handles name list,
#' still very preliminary
#'
#' @param x The iterated list item to be checked
#'
#' @return A formatted value to filter instruments with
#'
#' @noRd
InstrumentNameCheck <- function(x) {
  if (is.character(x)) {
    if (grepl("\\d", x) && grepl("[A-Za-z]", x)) {
      # message("Handle both")
      TheY <- x
    } else if (grepl("\\d", x)) {
      # message("Converting numeric")
      TheY <- as.numeric(x)
    } else if (grepl("[A-Za-z]", x)) {
      # message("Handle no letters")
      TheY <- x
    } else {
      stop("Instrument list item not character or numeric")
    }
  } else if (is.numeric(x)) {
    # message("Numeric")
    TheY <- x
  } else {
    stop("Instrument list item not character or numeric")
  }
  return(TheY)
}