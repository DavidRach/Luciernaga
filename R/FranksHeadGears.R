#' Internal FrankensteinsConfig,returns overall .xml config
#'
#' @param NumberDetectors Number detectors of alternate instrument
#'
#' @noRd
FranksHeadGears <- function(NumberDetectors) {
  InstrumentAndDetectors <- paste("Frankenstein", NumberDetectors, sep = " ")

  Metadata <- sprintf('<?xml version="1.0" encoding="utf-8"?>
<InstrumentConfiguration
  Name="%s"
  CytometerName="SpectralProt1"
  Institution="Cytekbio"
  SerialNumber="16888"
  BitsOfData="22"
  HasPlateLoader="True"
  HeaderSize="20">
  #PlaceLasersHere
</InstrumentConfiguration>', InstrumentAndDetectors)

  return(Metadata)
}