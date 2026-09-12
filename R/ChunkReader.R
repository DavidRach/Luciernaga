#' Internal for CytekQCFilePrep
#'
#' @importFrom utils read.csv
#'
#' @return An internal value
#'
#' @noRd
ChunkReader <- function(x){
  ReadChunks <- read.csv(x, check.names = FALSE)
}