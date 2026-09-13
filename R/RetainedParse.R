#' Internal for SC Unmix
#'
#' @param x Passed Retained Detector
#' @param data The Raw and Normed Data.frame
#'
#' @importFrom dplyr filter select pull
#' @importFrom tidyselect all_of
#' @importFrom rlang .data
#'
#' @return A vector of values to be passed to quantile
#' @noRd
RetainedParse <- function(x, data) {
  RetainedA <- paste0(x, "-A")
  MaxValues <- data %>%
    filter(.data[[x]] == 1) %>%
    select(all_of(RetainedA)) %>%
    pull(.) # TODO: `.` is explicit first-arg placeholder (magrittr-only)
}