
#' Internal for Simulated Data, coordinates pop level and adds up across fluorophore contributions
#'
#' @param x Iterated Population
#' @param ToAssemble Combined Population data
#' @param ScaledData The scaled signatures
#' @param distribution The data.frame with desired distribution parameters.
#'
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom dplyr pull
#' @importFrom purrr map
#' @importFrom purrr reduce
#'
#' @return The matrices of raw values for respective population
#'
#' @noRd
DataSimulation <- function(x, ToAssemble, ScaledData, distribution){
  IntData <- ToAssemble %>% filter(Pops %in% x) %>% select(-Pops)
  LocalNumber <- IntData %>% pull(Total)
  IntData <- IntData %>% select(-Total)
  NotZero <- IntData[, apply(IntData, 2, function(col) any(col != 0)), drop = FALSE]
  TheseFluorophores <- NotZero %>% colnames(.)

  if (length(TheseFluorophores) >= 1){
  #x <- TheseFluorophores[1]
  FluorMatrices <- map(.x=TheseFluorophores, .f=FluorophoreVariance,
   LocalNumber=LocalNumber, ScaledData=ScaledData, distribution=distribution)

  if (length(FluorMatrices) > 1){
    SummedContents <- reduce(FluorMatrices, `+`)
  } else {SummedContents <- FluorMatrices[[1]]}

  return(SummedContents)
  } else {message("Institute blank data protocol here")}
}