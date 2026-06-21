#' Small helper function, attempts to guess the correct max fluorophore for a
#' given instrument.
#' 
#' @param x The iterated in fluorophore name
#' @param NumberDetectors Value used by InstrumentReferences to identify an instrument
#' 
#' @importFrom dplyr filter arrange desc slice pull
#' 
#' @return If identified, the peak detector for a given fluor, else NULL
#' 
#' @export
#' 
#' @examples PeakDetector <- DetectorGuess(x="APC-Fire 810", NumberDetectors=64)
#' 
DetectorGuess <- function(x, NumberDetectors=64){

    TheKnownFluors <- Luciernaga:::InstrumentReferences(NumberDetectors=NumberDetectors)
    TheFluor <- TheKnownFluors |> filter(Fluorophore %in% x)

    if (nrow(TheFluor > 1)){
        TheDetector <- TheFluor |> arrange(desc(AdjustedY)) |> slice(1) |> pull(Detector)
    } else {TheDetector <- ""}

    return(TheDetector)
}