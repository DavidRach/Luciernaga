#' Internal for CytoSetCheck, generates ID
#' 
#' @param x Iterated in path to individual fcs files
#' 
#' @importFrom flowWorkspace load_cytoset_from_fcs
#' 
#' @return A concatinated string corresponding column name order
#' 
#' @noRd
CytoSetScreenInternal <- function(x){
  TheCytoSet <- load_cytoset_from_fcs(files=x)
  Values <- colnames(TheCytoSet[[1]])
  ID <- paste(Values, collapse = " ")
  ID <- data.frame(ID)
  return(ID)
}