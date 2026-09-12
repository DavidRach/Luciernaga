
#' Internal for Utility_Concatinate
#'
#' @param x Iterated metadata column name to edit on
#' @param data The concatenated data.frame to be edited
#' @param conversion The output of DoWeConvert
#'
#' @importFrom dplyr filter
#' @importFrom dplyr pull
#' @importFrom dplyr mutate
#' @importFrom dplyr across
#' @importFrom tidyselect all_of
#'
#' @noRd
ExecuteNumerics <- function(x, data, conversion){


  Numbers <- conversion  %>% dplyr::filter(Column %in% x) %>% pull(Numbers)
  Letters <- conversion  %>% dplyr::filter(Column %in% x) %>% pull(Letters)
  if (Numbers == TRUE && Letters==FALSE){
    data <- data %>% mutate(across(all_of(x), ~ as.numeric(.)))
  }

  if (Numbers == TRUE && Letters==TRUE){
    data <- data %>%
      mutate(across(all_of(x), ~ gsub("[A-Za-z]", "", .))) %>%
      mutate(across(all_of(x), ~ as.numeric(.)))
  }
  return(data)
}