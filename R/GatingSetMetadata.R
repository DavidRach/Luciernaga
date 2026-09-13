#' Internal Luciernaga_FolderGroupHeatmap, parses GatingSet metadata
#' 
#' @param x An iterated GatingSet object
#' 
#' @importFrom flowCore keyword
#' 
#' @return A one-row data.frame of GUID, TUBENAME, Date, and Count
#' 
#' @noRd
GatingSetMetadata <- function(x) {
  GUID <- keyword(x, "GUID")
  TUBENAME <- keyword(x, "TUBENAME")
  Date <- keyword(x, "$DATE")
  Count <- unname(nrow(x))[[1]]
  Data <- data.frame(GUID, TUBENAME, Date, Count)
  return(Data)
}