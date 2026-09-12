#' Internal for Stacked Reports
#'
#' @importFrom dplyr filter pull
#'
#' @noRd
ClusterAbundance <- function(x, data){
  #x <- TheClusters[1]

  Subset <- data |> dplyr::filter(Cluster %in% x)
  TheValues <- Subset |> pull(cutoff) |> unique()

  if(length(TheValues) == 1 && TheValues== TRUE){Value <- x
  } else if(length(TheValues) == 2 && any(TheValues==TRUE)){Value <- x
  } else {Value <- NULL}

 return(Value)
}