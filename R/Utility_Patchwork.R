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
#' @param NotListofList Internal default is True, don't touch.
#'
#' @importFrom purrr map
#' @importFrom patchwork wrap_plots
#' @importFrom grDevices dev.off
#' @importFrom grDevices pdf
#'
#' @return Some additional value to edit
#' @export
#'
#' @examples
#'
#' library(flowCore)
#' library(flowWorkspace)
#' library(openCyto)
#' library(data.table)
#'
#' File_Location <- system.file("extdata", package = "Luciernaga")
#' FCS_Files <- list.files(path = File_Location, pattern = ".fcs",
#'   full.names = TRUE)
#' UnstainedFCSFiles <- FCS_Files[grep("Unstained", FCS_Files)]
#' UnstainedCells <- UnstainedFCSFiles[-grep("Beads", UnstainedFCSFiles)]
#' MyCytoSet <- load_cytoset_from_fcs(UnstainedCells[1:3],
#'   truncate_max_range = FALSE,transformation = FALSE)
#' MyGatingSet <- GatingSet(MyCytoSet)
#' MyGates <- fread(file.path(path = File_Location, pattern = 'Gates.csv'))
#' MyGatingTemplate <- gatingTemplate(MyGates)
#' gt_gating(MyGatingTemplate, MyGatingSet)
#' removestrings <-  c("DR_", "Cells", ".fcs", "-", " ")
#' StorageLocation <- file.path("C:", "Users", "JohnDoe", "Desktop")
#'
#' AllSpecimens <- Utility_IterativeGating(x=MyGatingSet[1:3],
#'  sample.name = "GUID", removestrings = removestrings,
#'  subset = "nonDebris", gate = "lymphocytes", xValue = "FSC-A",
#'   yValue = "SSC-A", bins = 270)
#'
#' Patchwork <- Utility_Patchwork(AllSpecimens, "LymphocyteGates",
#'  outfolder=StorageLocation, thecolumns=2, therows=2,
#'   width = 7, height = 9, returntype="patchwork")
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
#'
#' @return An internal value
#'
#' @noRd
split_list <- function(input_list, chunk_size) {
  split(input_list, ceiling(seq_along(input_list) / chunk_size))
}

#' Wraps the sublist to hand off to patchwork
#' @importFrom patchwork wrap_plots
#'
#' @return An internal value
#'
#' @noRd
sublist_plots <- function(x, thecolumns, therows){
  p <- wrap_plots(x, ncol = thecolumns, nrow = therows, widths = 0.8, heights = 0.8)
}
