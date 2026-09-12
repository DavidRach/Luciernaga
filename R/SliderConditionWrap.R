
#' Internal for Wetlab_Decisions
#'
#' @param x The remaining conditions
#' @param y The remaining conditions cell aliquots
#' @param name The passed name
#' @param Date The passed date
#' @param TotalCells The passed TotalCells
#' @param RestConcentration The passed RestConcentration
#' @param FinalConcentration The passed FinalConcentration
#' @param TheConditions The Passed csv data
#'
#' @return An internal value
#'
#' @noRd
SliderConditionWrap <- function(x, y, name, Date, TotalCells, RestConcentration, FinalConcentration,
                                TheConditions){

  Condition <- x
  CellsPerTube <- y

  #if (!FinalConcentration==CellsPerTube){message("Cells are scarce resource for ", name)}

  FinalVolumeML <- (CellsPerTube*1)/FinalConcentration
  FinalVolumeML <- round(FinalVolumeML, 2)

  RestVolToAddML <- CellsPerTube/RestConcentration
  RestVolToAddML <- round(RestVolToAddML, 2)

  MediaVolToAddML <- FinalVolumeML-RestVolToAddML
  MediaVolToAddML <- round(MediaVolToAddML, 2)

  if (MediaVolToAddML < 0){SpinDown <- FALSE
  } else {SpinDown <- TRUE}

  TotalCells <- format(TotalCells, scientific=TRUE, digits=2)
  CellsPerTube <- format(CellsPerTube, scientific=TRUE, digits=2)
  FinalConcentration <- format(FinalConcentration, scientific=TRUE, digits=2)

  PreliminaryData <- cbind(name, Date, Condition, TotalCells, RestConcentration, FinalConcentration,
                           RestVolToAddML, MediaVolToAddML, CellsPerTube, FinalVolumeML)

  PreliminaryData <- data.frame(PreliminaryData, check.names=FALSE)
  return(PreliminaryData)
}