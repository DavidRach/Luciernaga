#' Internal for Luciernaga_Tree
#'
#' @param x An iterated path to a .csv to be read.
#'
#' @param utils read.csv
#' @param dplyr mutate
#' @param dplyr relocate
#'
#' @return An internal value
#'
#' @noRd
CSVRead <- function(x){

  name <- basename(x)
  internalstrings <- c("RelativeBrightness", ".csv")
  name <- Luciernaga::NameCleanUp(name, removestrings=internalstrings)
  Data <- read.csv(x, check.names=FALSE)
  Data <- Data %>% mutate(sample = name) %>% relocate(sample, .before=Cluster)
}