
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
FCS_Subset_Copy <- function(x, keywords, outpath, addon, metadataOverride){

if (!is.null(metadataOverride)){

    WhatTheHellMate <- do.call(
        rbind,
        lapply(keywords, function(x) {
            parts <- strsplit(sub("^~", "", x), "=")[[1]]
            data.frame(
            keyword = parts[1],
            value   = parts[2],
            stringsAsFactors = FALSE
            )}))
    
NotAsteriskArchive <- WhatTheHellMate[WhatTheHellMate$value != "*", ]
    
not_asterisk <- WhatTheHellMate$value != "*"
is_asterisk <- WhatTheHellMate$value == "*"

if (any(is_asterisk)) {
  resolved <- FlowKeywords(
    x        = x,
    keywords = WhatTheHellMate$keyword[is_asterisk],
    addon    = addon
  )
  WhatTheHellMate$value[is_asterisk] <- resolved
}
FileName <- paste(WhatTheHellMate$value, collapse = "_")
} else {FileName <- FlowKeywords(x=x, keywords=keywords, addon=addon)}

FileNameOut <- paste0(FileName, ".fcs")
if (is.null(outpath)) {outpath <- getwd()}
FinalRestingPlace <- file.path(outpath, FileNameOut)

new_fcs <- cytoframe_to_flowFrame(x)
    
if (!is.null(metadataOverride)){
    for (i in seq_len(nrow(NotAsteriskArchive))) {
    new_fcs@description[[NotAsteriskArchive$keyword[i]]] <-
    NotAsteriskArchive$value[i]
  }
}

new_fcs@description$GUID <- FileName
write.FCS(new_fcs, filename = FinalRestingPlace, delimiter="#")
}