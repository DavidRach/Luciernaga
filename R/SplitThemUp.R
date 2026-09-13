#' Internal for Chorizo, determines the Parameter row name order
#'
#' @param TheIntegers The iterated number of row names
#' @param RangeStart Default 10, because it matches
#' @param RangeSize Default 10, because it matches
#'
#' @importFrom purrr flatten
#'
#' @return A vector of dollar-p-number names to rearrange parameter data
#'
#' @noRd
SplitThemUp <- function(TheIntegers, RangeStart, RangeSize) {
  TheIntegers <- TheIntegers[TheIntegers >= RangeStart]
  RangesList <- split(TheIntegers,
                       (TheIntegers - RangeStart) %/% RangeSize + 1)
  RangesLength <- length(RangesList)
  InitialAppend <- 1:RangesLength
  Remainder <- 1:9
  These <- setdiff(Remainder, InitialAppend)

  for (i in seq_along(InitialAppend)) {
    RangesList[[i]] <- c(RangesList[[i]], InitialAppend[i])
  }

  RangesList[[RangesLength]] <- c(RangesList[[RangesLength]], These)

  RangesList <- flatten(RangesList)
  TheList <- unlist(RangesList)
  TheList <- paste0("$P", TheList)

  return(TheList)
}