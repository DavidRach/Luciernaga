#' Internal for Lucierna_Lists
#'
#' @param ListOfList A list of list
#' @param n Passed number of indices in the above
#'
#' @return An internal value
#'
#' @noRd
ItemSelect <- function(ListOfList, n) {
  result <- lapply(ListOfList, function(innerList) innerList[[n]])
  return(result)
}