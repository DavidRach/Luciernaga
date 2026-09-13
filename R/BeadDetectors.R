#' Internal for LuciernagaQC
#'
#' @param x The detector name to filter and summarize for
#' @param data The data.frame containing detector columns
#'
#' @importFrom dplyr filter select summarise
#' @importFrom stats median
#' @importFrom tidyselect all_of
#' @importFrom rlang .data
#'
#' @return An internal value
#'
#' @noRd
BeadDetectors <- function(x, data) {

  y <- paste0(x, "-A")
  FuckOff <- data |> filter(.data[[x]] == 1) |>
    select(all_of(y))
  TheDetector <- x
  colnames(FuckOff)[1] <- "Detector"
  TheMedian <- FuckOff |>
    summarise(TheMedian = median(Detector, na.rm = TRUE))
  Return <- cbind(TheDetector, TheMedian)
}