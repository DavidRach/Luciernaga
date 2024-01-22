#' Originial Subtract Autofluorescence from Median Single Colors
#'
#' @param i A row of the data.frame object
#' @param Data A data.frame object
#' @importFrom dplyr filter
#'
#' @return NULL
#' @export
#'
#'
#' @examples NULL
Utility_Subtract <- function(i, Data){
  #Generalization Note: This one may be a lot of effort generalize, most of it is specific.

  Z <- Data[i,]
  Alpha <- Data %>% filter(Experiment == Z$Experiment) %>% filter(Type == Z$Type) %>% filter(Negative == Z$Negative) %>% filter(str_detect(Ligand, "Unstained"))
  Subtracted <- Z[,6:69] - Alpha[,6:69]
  Updated <- cbind(Z[,1:5], Subtracted)
  Updated %>% mutate(Processed = "Yes") %>% relocate(Processed, .after = Negative)
}
