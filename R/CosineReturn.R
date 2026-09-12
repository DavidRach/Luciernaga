#' Internal for InstrumentComparison, derrives the cosine data
#' for the fluors and instruments
#'
#' @param x The instrument being filtered for
#' @param MainFluorophore The fluorophore the rest are being compared to
#' @param data The data.frame containing required data
#'
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom tidyselect where
#'
#' @return A single column renamed for the instrument with
#'  comparisons as rows
#'
#' @noRd
CosineReturn <- function(x, MainFluorophore, data){
  Subset <- data |> dplyr::filter(Instrument %in% x)
  Poised <- Subset |> select(-Instrument)
  Poised2 <- Poised %>% select(where(~ !all(is.na(.))))
  CosineReturn <- Luciernaga_Cosine(data=Poised2,
   returntype="data", rearrange = FALSE)
  CosineReturn <- round(CosineReturn, 2)
  Data <- data.frame(CosineReturn, check.names=FALSE)
  TheValues <- Data |> select(MainFluorophore)
  colnames(TheValues)[1] <- x
  return(TheValues)
}