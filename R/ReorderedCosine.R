#' Internal for LuciernagaReport
#'
#' @importFrom stats as.dist hclust
#'
#' @return An internal value
#'
#' @noRd
ReorderedCosine <- function(CosineMatrix){
  Day <- as.dist((1-CosineMatrix)/2)
  Night <- hclust(Day)
  Twilight <- CosineMatrix[Night$order, Night$order]
}