#' Downsample a gs based on subset
#'
#' @param x A gating set object
#' @param sample.name Keyword specifying sample name
#' @param removestrings Value to be removed from sample name
#' @param subsets The gating hierarchy subset you want to include
#' @param subsample Total number of events to sub-sample from each specimen
#' @param inverse.transform Whether to reverse the GatingSet Transform on the data,
#' default is set to FALSE.
#' @param internal Whether to return as a data.frame (vs. a flow frame or .fcs file)
#' @param export Default is set to FALSE, when TRUE (and internal = FALSE) returns a
#' .fcs file to outpath.
#' @param outpath When export is true, the file.path to where you want the .fcs file
#' stored.
#'
#' @importFrom flowWorkspace keyword
#' @importFrom flowWorkspace gs_pop_get_data
#' @importFrom flowCore exprs
#' @importFrom dplyr slice_sample
#' @importFrom flowCore parameters
#' @importFrom flowCore write.FCS
#' @importFrom dplyr mutate
#' @importFrom methods new
#'
#' @return Either a data.frame, a flow.frame or an .fcs file depending on your
#' selected options
#' @export
#'
#' @examples NULL


Utility_Downsample <- function(x, sample.name, removestrings,
                                  subsets, subsample=NULL, inverse.transform = FALSE,
                                  internal = FALSE, export = FALSE, outpath = NULL){

  name <- keyword(x, sample.name)
  alternatename <- NameCleanUp(name, removestrings)

  #Retrieving the exprs data for my subset population of interest
  ff <- gs_pop_get_data(x, subsets, inverse.transform = inverse.transform)
  df <- exprs(ff[[1]])
  DF <- as.data.frame(df, check.names = FALSE)

  # If down-sampling is specified
  if (!is.null(subsample)) {DF <- slice_sample(DF, n = subsample,
                                              replace = FALSE)
  } else {DF <- DF}

  if (internal == FALSE){
    FCSSubset <- as.matrix(DF)
    FlowFrameTest <- ff[[1, returnType = "flowFrame"]]
    original_p <- parameters(FlowFrameTest)
    original_d <- keyword(FlowFrameTest)
    new_fcs <- new("flowFrame", exprs=FCSSubset, parameters=original_p,
                   description=original_d)

    if (export == TRUE) {
      if (!is.null(outpath)) {TheFileName <- paste0(alternatename, "_downsampled.fcs")
                              fileSpot <- file.path(outpath, TheFileName)}
                              write.FCS(new_fcs, filename = fileSpot, delimiter="#")
    } else {return(new_fcs)}

  } else { message(alternatename, " has been processed")
           DF <- DF %>% mutate(specimen = alternatename)
           return(DF)}
}
