#' Small Internal Function
#' 
#' @param data Something
#' @param x Something
#' @param type Something
#' 
#' @importFrom dplyr select
#' @importFrom tidyselect all_of
#' @importFrom dplyr slice
#' @importFrom dplyr pull
#' @importFrom stats na.omit
#' @return Some value
#' 
#' @noRd
CurrentStatus <- function(data, x, type){
  Status <- data |> select(all_of(x)) |> na.omit() |>
    slice(1)|> pull()
  return(Status)
}