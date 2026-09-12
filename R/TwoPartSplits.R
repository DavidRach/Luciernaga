#' Helper function for QC PDF conversions. Tidyes two chunk items into single data.frame cell
#' 
#' @param x A string line read in from the pdf
#' 
#' @noRd
TwoPartSplits <- function(x){
    x <- gsub("mum:", "mum:  ", x)
    x <- gsub("High", "    High", x)
    parts <- strsplit(trimws(x), "\\s{2,}")[[1]]
    keys <- trimws(gsub(":$", "", parts[seq(1, length(parts), by = 2)]))
    vals <- trimws(parts[seq(2, length(parts), by = 2)])
    if (length(keys) > length(vals)) {
    keys <- keys[seq_len(length(vals))]}
    data <- as.data.frame(as.list(setNames(vals, keys)), check.names=FALSE)
    return(data)
}