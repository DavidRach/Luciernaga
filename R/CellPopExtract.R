#' Extracts out a subset from a GatingSet, and returns as their respective .fcs files
#'  to a target folder. Additional option to downsample to a desired number. 
#' 
#' @param x Either a file path to a FlowJo workspace, or a GatingSet object. 
#' @param path When using a FlowJo workspace for argument x, path to the folder containing
#' the respective .fcs files contained within that workspace
#' @param keywords Default "GROUPNAME", can handle up to 3, ex. c("$PROJ", "GROUPNAME", "TUBENAME")
#' @param subset The population node in the GatingSet that you want to extract as their own .fcs files
#' @param outpath The file.path to the storage location
#' @param addon Default NULL, appends to the end of the filename to distinguish from original fcs file.
#' @param downsample Default NULL, provide a number to downsample to that number, alternatively use
#' 0.1 for a proportion of total cells
#' 
#' @importFrom CytoML open_flowjo_xml flowjo_to_gatingset
#' @importFrom flowWorkspace gs_pop_get_data
#' @importFrom purrr walk
#' 
#' @return Returns .fcs files for the designated cell population to the designated folder. 
#' 
#' @export
#' 
#' @examples A <- 2+2
CellPopExtract <- function(x, path, keywords="GROUPNAME", subset, outpath, addon=NULL, 
downsample=NULL){

if (inherits(x, "GatingSet")){gs <- x
} else {
ws <- CytoML::open_flowjo_xml(x)
gs <- CytoML::flowjo_to_gatingset(ws, name=1, path = path, keywords = keywords)
}

if (!is.null(downsample)){
    purrr::walk(.x=gs, .f=Utility_Downsample, sample.name = keywords,
    subsets = subset, subsample = downsample, internal = FALSE, export = TRUE,
    inverse.transform=TRUE, outpath = outpath, addon=addon)
} else {
    SubsetSolo <- flowWorkspace::gs_pop_get_data(gs, subset, inverse.transform = TRUE)
    purrr::walk(.x=SubsetSolo, .f=FCS_Subset_Copy, keywords=keywords, 
        outpath=outpath, addon=addon)
}
}


#' Internal for CellPopExtract
#' 
#' @param x Either a file path to a FlowJo workspace, or a GatingSet object. 
#' @param keywords Default "GROUPNAME", can handle up to 3, ex. c("$PROJ", "GROUPNAME", "TUBENAME")
#' @param outpath The file.path to the storage location
#' @param addon Default NULL, appends to the end of the filename to distinguish from original fcs file.
#' 
#' @importFrom flowCore keyword write.FCS
#' @importFrom flowWorkspace cytoframe_to_flowFrame
#' 
#' @return Writes the .fcs file with altered naming to the designated outpath
#' 
#' @noRd
FCS_Subset_Copy <- function(x, keywords, outpath, addon){

FileName <- FlowKeywords(x=x, keywords=keywords, addon=addon)

FileNameOut <- paste0(FileName, ".fcs")
if (is.null(outpath)) {outpath <- getwd()}
FinalRestingPlace <- file.path(outpath, FileNameOut)

new_fcs <- cytoframe_to_flowFrame(x)
new_fcs@description$GUID <- FileName
write.FCS(new_fcs, filename = FinalRestingPlace, delimiter="#")
}

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