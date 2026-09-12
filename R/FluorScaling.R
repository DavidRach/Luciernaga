#' Internal for SimulatedData, adjust reference signatures to desired MFI targets
#'
#' @param x Iterated in Fluorophore
#' @param data Iterated in Fluorophore Signature Data
#' @param targets The data.frame containing Fluorophore MFI targets
#'
#' @importFrom dplyr filter
#' @importFrom dplyr pull
#' @importFrom dplyr mutate
#'
#' @return MFI-scaled Fluorophore Signature
#'
#' @noRd
FluorScaling <- function(x, data, targets){
  IntData <- data %>% filter(Fluorophore %in% x)
  Multiplier <- targets %>% filter(Fluorophores %in% x) %>% pull(MFI)
  IntData <- IntData %>% mutate(AdjustedY = AdjustedY*Multiplier)
  return(IntData)
}