#' Splits available plots into sublist
#'
#' @return An internal value
#'
#' @noRd
split_list <- function(input_list, chunk_size) {
  split(input_list, ceiling(seq_along(input_list) / chunk_size))
}