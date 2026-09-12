#' Helper function for QC PDF conversions. Combines the bread ends,
#'  then separates out the middle as it's own separate data.frame,
#'  returning both as components of a list. 
#' 
#' @param data The lines from the pdf for the respective section
#' 
#' @noRd
BridgeSplits <- function(data){
    Total <- length(data)
    parts <- strsplit(trimws(data), "\\s{2,}")[[1]]
    Second <- paste0(parts[1], ": ", parts[3])
    Bridge <- as.data.frame(t(setNames(strsplit(Second, ":\\s*")[[1]][2],
                               strsplit(Second, ":\\s*")[[1]][1])),
                                check.names=FALSE)
    TheParts <- strsplit(trimws(data), "\\s{2,}")[2:6]
    TheVector <- c(parts[2], unlist(TheParts))

    BridgeData <- do.call(rbind, lapply(TheVector, function(s) {
    parts <- strsplit(s, ":\\s*")[[1]]
    data.frame(
        Name  = parts[1],
        Value = as.numeric(parts[2]),
        stringsAsFactors = FALSE
    )
    }))

    ReturnVals <- list(Bridge, BridgeData)
    return(ReturnVals)
}