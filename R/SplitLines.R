
#' Internal for PDF Processing, handles Chorus 7 column QC params
#' 
#' @param x A read-in line from the pdf file
#' 
#' @noRd
SplitLines <- function(x){
    parts <- strsplit(trimws(x), "\\s{2,}")[[1]]
    if (!(length(parts) %in% c(3, 7))){parts <- NULL}
    return(parts)
}