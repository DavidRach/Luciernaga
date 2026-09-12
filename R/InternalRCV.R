#' Internal for RCV from FCS
#' 
#' @param x The column name for an individual marker
#' @param data The data.frame of exprs values
#' 
#' @importFrom dplyr select
#' @importFrom stats mad
#' @importFrom stats median
#' 
#' @return An RCV for the individual marker
#' 
#' @noRd
InternalRCV <- function(x, data){
  Name <- x
  TheCol <- data |> select(x)
  TheCol <- TheCol[!is.na(TheCol)]
  if(length(TheCol) == 0) return(NULL)
  RCV <- mad(TheCol) / median(TheCol)
  RCV <- data.frame(RCV)
  colnames(RCV) <- Name
  return(RCV)
}