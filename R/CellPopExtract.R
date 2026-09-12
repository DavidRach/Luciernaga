#' Extracts out a subset from a GatingSet, and returns as their
#'  respective .fcs files to a target folder. Additional option to
#'  downsample to a desired number. 
#' 
#' @param x Either a file path to a FlowJo workspace, or a
#'  GatingSet object. 
#' @param path When using a FlowJo workspace for argument x,
#'  path to the folder containing
#' the respective .fcs files contained within that workspace
#' @param keywords Default "GROUPNAME", can handle up to 3,
#'  ex. c("$PROJ", "GROUPNAME", "TUBENAME")
#' @param subset The population node in the GatingSet that
#'  you want to extract as their own .fcs files
#' @param outpath The file.path to the storage location
#' @param addon Default NULL, appends to the end of the filename
#'  to distinguish from original fcs file.
#' @param downsample Default NULL, provide a number to downsample
#'  to that number, alternatively use 0.1 for a proportion of total
#'  cells
#' @param metadataOverride Default NULL
#' 
#' @importFrom CytoML open_flowjo_xml flowjo_to_gatingset
#' @importFrom flowWorkspace gs_pop_get_data
#' @importFrom purrr walk
#' 
#' @return Returns .fcs files for the designated cell population to
#'  the designated folder. 
#' 
#' @export
#' 
#' @examples A <- 2+2
CellPopExtract <- function(x, path, keywords="GROUPNAME", subset,
 outpath, addon=NULL, downsample=NULL, metadataOverride=NULL){

if (inherits(x, "GatingSet")){gs <- x
} else {
ws <- CytoML::open_flowjo_xml(x)
gs <- CytoML::flowjo_to_gatingset(ws, name=1, path = path,
     keywords = keywords)
}

if (!is.null(downsample)){
    purrr::walk(.x=gs, .f=Utility_Downsample,
         sample.name = keywords,
    subsets = subset, subsample = downsample, internal = FALSE,
     export = TRUE,
    inverse.transform=TRUE, outpath = outpath, addon=addon)
} else {
    SubsetSolo <- flowWorkspace::gs_pop_get_data(gs, subset,
         inverse.transform = TRUE)
    purrr::walk(.x=SubsetSolo, .f=FCS_Subset_Copy,
         keywords=keywords,  outpath=outpath, addon=addon,
         metadataOverride=metadataOverride)
}
}