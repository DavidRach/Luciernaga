#' Internal for Stacked Reports
#'
#' @param x The cluster identity to filter for
#' @param data The data.frame containing Cluster and cutoff columns
#'
#' @importFrom dplyr filter pull
#'
#' @return The cluster identity if its cutoff qualifies, else NULL
#'
#' @noRd
ClusterAbundance <- function(x, data) {

  #x <- TheClusters[1]

  Subset <- data |> filter(Cluster %in% x)
  TheValues <- Subset |> pull(cutoff) |> unique()

  if (length(TheValues) == 1 && TheValues== TRUE) {
    Value <- x
  } else if (length(TheValues) == 2 && any(TheValues==TRUE)) {
    Value <- x
  } else {
    Value <- NULL
  }

  return(Value)
}