#' Internal for Chorizo, rearranges dollar-p-number according to bizarre order
#'
#' @param x The parameter data.frame
#'
#' @return The rearranged parameter data.frame according to wacky order.
#' @noRd
RowNameReorder <- function(x) {
  TheRowNames <- rownames(x)
  LastElement <- TheRowNames[length(TheRowNames)]
  LastElement <- gsub("$P", "", fixed = TRUE, LastElement)
  LastNumber <- as.numeric(LastElement)

  TheIntegers <- 1:LastNumber

  TheOrder <- SplitThemUp(TheIntegers, RangeStart = 10, RangeSize = 10)

  Rearranged <- x[TheOrder, ]

  return(Rearranged)
}