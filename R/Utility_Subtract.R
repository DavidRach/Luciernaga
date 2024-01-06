#' Originial Subtract Autofluorescence from Median Single Colors
#'
#' @param i A row of the data.frame object
#' @param Data A data.frame object
#'
#' @return NULL
#' @export NULL
#'
#' @examples NULL
Utility_Subtract <- function(i, Data){
  Z <- Data[i,]
  Alpha <- Data %>% dplyr::filter(Experiment == Z$Experiment) %>% dplyr::filter(Type == Z$Type) %>% dplyr::filter(Negative == Z$Negative) %>% dplyr::filter(str_detect(Ligand, "Unstained"))
  Subtracted <- Z[,6:69] - Alpha[,6:69]
  Updated <- cbind(Z[,1:5], Subtracted)
  Updated %>% mutate(Processed = "Yes") %>% relocate(Processed, .after = Negative)
}
