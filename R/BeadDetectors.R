
#' Internal for LuciernagaQC
#'
#' @importFrom dplyr filter select summarise
#' @importFrom stats median
#' @importFrom tidyselect all_of
#'
#' @return An internal value
#'
#' @noRd
BeadDetectors <- function(x, data){

  y <- paste0(x, "-A")
  FuckOff <- data %>% dplyr::filter(.data[[x]] == 1) |>
    select(all_of(y))
  TheDetector <- x
  colnames(FuckOff)[1] <- "Detector"
  TheMedian <- FuckOff |>
    summarise(TheMedian = median(Detector, na.rm = TRUE))
  Return <- cbind(TheDetector, TheMedian)
}