#' Internal for Utility_Concatinate
#'
#' @param x The iterated metadata column being acted on.
#' @param data The concatinated data.frame
#' @param dictionary The reference table data.frame with the new numeric factors to swap out
#' the old character values
#'
#' @importFrom rlang !!
#' @importFrom dplyr select
#' @importFrom rlang sym
#' @importFrom rlang :=
#' @importFrom dplyr mutate
#' @importFrom dplyr recode
#' @importFrom stats setNames
#'
#' @noRd
ExecuteSwap <- function(x, data, dictionary){
  newname <- paste0("New_", x)
  data %>% select(!!sym(x))

  data <- data %>% mutate(!!sym(x) := recode(!!sym(x), !!!setNames(dictionary[[newname]],
                                                                   dictionary[[x]])))
  return(data)
}