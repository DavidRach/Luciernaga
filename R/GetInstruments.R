#' Internal Helper for InstrumentComparison
#' Returns codes for all instruments
#' 
#' @return List of Instrument Codes
#' 
#' @noRd
GetInstruments <- function() {
  Available <- c(64, 54, 48, 38, 30, 22, 14,
                 184, "182_DUV", 147, 112, 86,
                 "52_7L", 55, 51, "52_6L",
                 78, "48_A", "48_S")
  return(Available)
}