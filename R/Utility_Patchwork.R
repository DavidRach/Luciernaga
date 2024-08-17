#' Generate the plots into a single pdf
#'
#' @param x A list of plots
#' @param filename Name to save the .pdf as
#' @param outfolder Location to save the .pdf
#' @param thecolumns The number of columns per page
#' @param therows The number of rows per page
#' @param width Desired page width
#' @param height Desired page height
#' @param returntype Whether to return "pdf" (to desired location) or "patchwork" (to R)
#' @param NotListOfList Internal default is True, don't touch.
#'
#' @importFrom purrr map
#' @importFrom patchwork wrap_plots
#' @importFrom grDevices dev.off
#' @importFrom grDevices pdf
#'
#' @return Some additional value to edit
#' @export
#'
#' @examples NULL
#'

Utility_Patchwork <- function(x, filename, outfolder, thecolumns=2, therows=3,
  width = 7, height = 9, returntype="pdf", NotListofList = TRUE){

  if (NotListofList == TRUE){
  theList <- x
  theList <- Filter(Negate(is.null), theList)
  theListLength <- length(theList)

  theoreticalitems <- therows*thecolumns

  sublists <- split_list(theList, theoreticalitems)
  #length(sublists)
  } else{sublists <- x}

  if (returntype == "pdf"){

    MergedName <- paste(outfolder, filename, sep = "/")
    pdf(file = paste(MergedName, ".pdf", sep = "", collapse = NULL), width = width,
        height = height)
    p <- map(sublists, .f=sublist_plots, thecolumns=thecolumns, therows=therows)
    print(p)
    dev.off()

  } else if (returntype == "patchwork"){
      p <- map(sublists, .f=sublist_plots, thecolumns=thecolumns, therows=therows)
      return(p)}
}

#' Splits available plots into sublist
#' @noRd
split_list <- function(input_list, chunk_size) {
  split(input_list, ceiling(seq_along(input_list) / chunk_size))
}

#' Wraps the sublist to hand off to patchwork
#' @importFrom patchwork wrap_plots
#' @noRd
sublist_plots <- function(x, thecolumns, therows){
  p <- wrap_plots(x, ncol = thecolumns, nrow = therows, widths = 0.8, heights = 0.8)
}
