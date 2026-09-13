#' Internal for LuciernagaReport
#'
#' @param CosineMatrix TBD
#'
#' @importFrom stats as.dist hclust
#'
#' @return A reordered cosine matrix, clustered via hierarchical clustering
#'
#' @noRd
ReorderedCosine <- function(CosineMatrix) {
  Day <- as.dist((1 - CosineMatrix) / 2)
  Night <- hclust(Day)
  Twilight <- CosineMatrix[Night$order, Night$order]
}