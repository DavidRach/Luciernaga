
#' Internal for Utility_Concatinate
#'
#' @param x The iterated metadata column name
#' @param data The concatinated data.frame
#' @param conversion The reference dictionary data.frame table
#'
#' @importFrom dplyr filter
#' @importFrom dplyr pull
#' @importFrom dplyr select
#' @importFrom tidyselect all_of
#' @importFrom dplyr arrange
#' @importFrom dplyr desc
#' @importFrom dplyr mutate
#' @importFrom rlang !!
#' @importFrom rlang :=
#'
#' @noRd
ExecuteCharacters <- function(x, data, conversion){

  Numbers <- conversion  %>% dplyr::filter(Column %in% x) %>% pull(Numbers)
  Letters <- conversion  %>% dplyr::filter(Column %in% x) %>% pull(Letters)
  if (Numbers == FALSE && Letters==TRUE){
    newName <- paste0("New_", x)
    Internal <- data %>% select(all_of(x))
    SpecimenNames <- data.frame(table(Internal))
    #colnames(SpecimenNames)[[1]] <- "specimen"
    TheCount <- paste0("Count_", x)
    colnames(SpecimenNames)[[2]] <- TheCount
    SpecimenNames <- SpecimenNames %>% arrange(desc(TheCount))
    SpecimenNames <- SpecimenNames %>% mutate(!!newName := as.numeric(factor(x)))
    return(SpecimenNames)
  } else {SpecimenNames <- NULL}
}