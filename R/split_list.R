#' Splits available plots into sublist
#'
#' @param input_list TBD
#' @param chunk_size TBD
#'
#' @return An internal value
#'
#' @noRd
split_list <- function(input_list, chunk_size) {
  split(input_list, ceiling(seq_along(input_list) / chunk_size))
}