#' Concatenate a gs based on subset and optional subsample
#'
#' @param gs A gating set object
#' @param sample.name Keyword specifying sample name
#' @param removestrings Value to be removed from sample name
#' @param subsets The gating hierarchy subset you want to include
#' @param subsample Default is set to NULL. If down-sample is desired, Total number of events to subsample from each specimen
#' @param ReturnType Whether to return as a "data.frame", a "flow.frame", or export to an "fcs" file to the outpath location
#' @param newName File Name for the Concatenate File
#' @param outpath Location to store the concatenated file
#' @param export Whether to export as a .fcs file.
#' @param inverse.transform Whether to reverse the GatingSet Transform on the data, default is set to FALSE.
#'
#' @return A concatenated .fcs file to new location
#' @export
#'
#' @examples NULL

Utility_Concatinate <- function(gs, sample.name, removestrings, subsets, subsample=NULL, ReturnType, newName, outpath=NULL, export=FALSE,
                                inverse.transform = FALSE) {

gs <- gs
sample.name <- sample.name
removestrings <- removestrings
subsample <- subsample

if (!is.null(subsample)){
  if (!is.numeric(subsample)) {stop("Subsample argument is not numeric, enter just a number, no quotes needed")}
}

ConcatenatedFile <- map(gs, .f=Utility_Downsample, sample.name=sample.name, removestrings=removestrings,
                        subsets=subsets, subsample=subsample, internal = TRUE, export=FALSE) %>% bind_rows

if (ReturnType == "data.frame"){return(ConcatenatedFile)
}

ff <- gs_pop_get_data(gs, subsets, inverse.transform = inverse.transform)
ConcatenatedExtra <- ConcatenatedFile %>% select(specimen)
ConcatenatedMain <- ConcatenatedFile %>% select(-specimen)

new_fcs <- Utility_ColAppend(ff=ff, DF=ConcatenatedMain, columnframe=ConcatenatedExtra)

# Renaming of the .fcs file in the internal parameters would occur here.

if (ReturnType == "flow.frame") {return(new_FCS_2)}

if (ReturnTye  == "fcs"){
  if (is.null(outpath)) {outpath <- getwd()}
  TheFileName <- paste0(newName, ".fcs")
  fileSpot <- file.path(outpath, TheFileName)
  write.FCS(new_FCS_2, filename = fileSpot, delimiter="#")
}

}

