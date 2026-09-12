
#' Internal for Wetlab_Decisions
#'
#' @param x The passed conditions
#' @param name The passed name
#' @param Date The passed date
#' @param TotalCells The passed total cells
#' @param RestConcentration The passed rest concentration
#' @param FinalConcentration The passed final concentration
#' @param CellsPerTube The passed cells per tube
#' @param TheConditions The passed CSV data.
#'
#' @return An internal value
#'
#' @noRd
ConditionWrap <- function(x, name, Date, TotalCells,
  RestConcentration, FinalConcentration, CellsPerTube,
  TheConditions){

  Condition <- x

  FinalVolumeML <- (CellsPerTube*1)/FinalConcentration
  FinalVolumeML <- round(FinalVolumeML, 2)

  RestVolToAddML <- FinalConcentration/RestConcentration
  RestVolToAddML <- round(RestVolToAddML, 2)

  MediaVolToAddML <- FinalVolumeML-RestVolToAddML
  MediaVolToAddML <- round(MediaVolToAddML, 2)

  if (MediaVolToAddML < 0){SpinDown <- FALSE
  } else {SpinDown <- TRUE}

  FinalConcentration <- format(FinalConcentration,
     scientific=TRUE, digits=2)
  CellsPerTube <- format(CellsPerTube,
     scientific=TRUE, digits=2)

  PreliminaryData <- cbind(name, Date, Condition, TotalCells,
    RestConcentration, FinalConcentration, RestVolToAddML,
    MediaVolToAddML, CellsPerTube, FinalVolumeML)

  PreliminaryData <- data.frame(PreliminaryData, check.names=FALSE)
  return(PreliminaryData)
}