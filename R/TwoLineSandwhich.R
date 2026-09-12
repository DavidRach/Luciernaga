#' Helper function for QC PDF conversions. Combines the bread ends,
#'  then handles the middle, then tackles the middle continuation 
#'  on the subsequent row.
#' 
#' @param x A string line index from the pdf
#' @param data The lines from the pdf for filtering on
#' 
#' @importFrom dplyr bind_cols
#' 
#' @noRd
TwoLineSandwhich <- function(x, data){
    FirstLine <- data[x]
    SecondLine <- data[x+1]

    parts <- strsplit(trimws(FirstLine), "\\s{2,}")[[1]]
    First <- paste0(parts[1], " ", parts[2])
    Second <- paste0(parts[1], ": ", parts[3])
    parts <- c(First, Second)
    FirstData <- as.data.frame(
        setNames(
            list(
            trimws(sub(".*:","", parts[1])),
            trimws(sub(".*:","", parts[2]))
            ),
            trimws(sub(":.*","", parts))
        ), stringsAsFactors = FALSE, check.names=FALSE)

    SecondLine <- gsub("mum:", "mum:  ", SecondLine)
    SecondLine <- gsub("ity:", "ity:  ", SecondLine)
    SecondLine <- gsub("wer:", "wer:  ", SecondLine)
    parts <- strsplit(trimws(SecondLine), "\\s{2,}")[[1]]
    keys <- trimws(gsub(":$", "", parts[seq(1, length(parts), by = 2)]))
    vals <- trimws(parts[seq(2, length(parts), by = 2)])
    SecondData <- as.data.frame(as.list(setNames(vals, keys)), check.names=FALSE)
    colnames(SecondData) <- paste0(colnames(FirstData[2]), " ", colnames(SecondData))
    TheData <- bind_cols(FirstData, SecondData)
    return(TheData)
}