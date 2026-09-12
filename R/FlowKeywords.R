
#' An util function, handles the sample.name/keyword piping, handling generating
#' an appended name out of flowCore keywords
#' 
#' @param x The iterated GatingSet object
#' @param keywords The .fcs file keyword, or a c("A", "B", "C") vector of keywords to extract
#' @param addon A character string value to append after the x argument, but before the .fcs
#' 
#' @importFrom flowCore keyword
#' 
#' @return A character string resulting from the provided keywords for the respective specimen. 
#' 
#' @noRd
FlowKeywords <- function(x, keywords, addon){
    if (length(keywords) > 0){
        First <- keyword(x, keywords[1])

    if (length(keywords) > 1) {
        Second <- keyword(x, keywords[2])
    
    if (length(keywords) > 2) {
        Third <- keyword(x, keywords[3])
   
    if (length(keywords) > 3) {
        stop("Please choose only three keywords, thank you!")
    } else {Nomenclature <- paste(First, Second, Third, sep="_")} 
    } else {Nomenclature <- paste(First, Second, sep="_")}    
    } else {Nomenclature <- First}
    } else {stop("No keywords provided")}

    if (!is.null(addon)){FileName <- paste(Nomenclature, addon, sep="_")
        } else {FileName <- Nomenclature}
    
    return(FileName)
}