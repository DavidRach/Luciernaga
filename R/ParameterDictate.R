#' Internal for Chorizo, iterates out individual detector parameters
#'
#' @param x Iterated In Dollar-P-Number for filtering
#' @param data The reordered parameter data
#'
#' @importFrom dplyr pull
#' @importFrom stringr str_detect
#'
#' @return The individual detectors list of parameters for description file
#' @noRd
ParameterDictate <- function(x, data) {
  Subset <- data[x, ]

  Bits <- paste0(x, "B")
  Ehh <- paste0(x, "E")
  Name <- paste0(x, "N")
  Range <- paste0(x, "R")
  Type <- paste0(x, "TYPE")
  Volt <- paste0(x, "V")

  BitsVal <- "32"
  EhhVal <- "0,0"
  NameVal <- Subset |> pull(name)
  RangeVal <- Subset |> pull(range)
  RangeVal <- as.character(RangeVal)

  if (str_detect(NameVal, "FSC")) {
    TypeVal <- "Forward_Scatter"
  } else if (str_detect(NameVal, "SSC")) {
    TypeVal <- "Side_Scatter"
  } else if (str_detect(NameVal, "Time")) {
    TypeVal <- "Time"
  } else {
    TypeVal <- "Raw_Fluorescence"
  }

  if (!str_detect(NameVal, "Time")) {
    VoltVal <- "307"
    DetectorList <- list(
      Bits = BitsVal,
      Ehh = EhhVal,
      Name = NameVal,
      Range = RangeVal,
      Type = TypeVal,
      Volt = VoltVal
    )
    names(DetectorList) <- c(Bits, Ehh, Name, Range, Type, Volt)
  } else {
    DetectorList <- list(
      Bits = BitsVal,
      Ehh = EhhVal,
      Name = NameVal,
      Range = RangeVal,
      Type = TypeVal
    )
    names(DetectorList) <- c(Bits, Ehh, Name, Range, Type)
  }

  return(DetectorList)
}