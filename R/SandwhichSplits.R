#' Helper function for QC PDF conversions. Combines the bread ends, then handles the middle.
#' 
#' @param x A string line read in from the pdf
#' 
#' @importFrom stats setNames
#' 
#' @return A one-row data.frame with two named columns parsed from the line
#' 
#' @noRd
SandwhichSplits <- function(x) {
  parts <- strsplit(trimws(x), "\\s{2,}")[[1]]
  First <- paste0(parts[1], " ", parts[2])
  Second <- paste0(parts[1], ": ", parts[3])
  parts <- c(First, Second)
  data <- as.data.frame(
    setNames(
      list(
        trimws(sub(".*:", "", parts[1])),
        trimws(sub(".*:", "", parts[2]))
      ),
      trimws(sub(":.*", "", parts))
    ), stringsAsFactors = FALSE, check.names = FALSE)
  return(data)
}