#' Internal for Wetlab_Decisions
#'
#' @param x The condition insufficient
#' @param name The passed name
#' @param Date The passed date
#' @param TotalCells The passed Total Cells
#' @param RestConcentration The passed rest concentration
#'
#' @importFrom stats setNames
#'
#' @return An internal value
#'
#' @noRd
TheInsufficient <- function(x, name, Date, TotalCells, RestConcentration) {

  Condition <- x
  FinalConcentration <- NA
  RestVolToAddML <- NA
  MediaVolToAddML <- NA
  CellsPerTube <- NA
  FinalVolumeML <- NA

  PreliminaryData <- cbind(
    name, Date, Condition, TotalCells, RestConcentration, FinalConcentration,
    RestVolToAddML, MediaVolToAddML, CellsPerTube, FinalVolumeML
  )

  PreliminaryData <- data.frame(PreliminaryData, check.names = FALSE)
  return(PreliminaryData)
}