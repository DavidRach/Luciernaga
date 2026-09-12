#' Parses the two-page Amnis FlowSight QC report, returning tidy dataframe.
#'
#' @param x A file.path to the respective QC report pdf
#' 
#' @importFrom pdftools pdf_text
#' @importFrom dplyr bind_cols
#' 
#' @export
#' 
#' @examples A <- 2+2
QC_FlowSightPDF <- function(x){
  text <- pdftools::pdf_text(x)
  NumberPages <- length(text)
  PageOne <- FlowSightPageOne(x=text[1])
  PageTwo <- FlowSightPageTwo(x=text[2])
  Completed <- bind_cols(PageOne, PageTwo)
  return(Completed)
}