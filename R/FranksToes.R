#' Internal FrankensteinsConfig, generates detector outputs
#' 
#' @param x The iterated Detector
#' @param detector The Detector data.frame
#' 
#' @importFrom dplyr filter mutate pull
#' 
#' @noRd
FranksToes <- function(x, detector){
  Target <- detector |> filter(Detectors %in% x)

  NameValue <- Target |> mutate(Name = Detectors) |> pull(Name)
  NameValue <- sub(".*-", "", NameValue)
  TheWave <- sub("/.*", "", NameValue)
  TheBand <- sub(".*/", "", NameValue)

  Name <- Target |> pull(Detectors)
  Number <- Target |> pull(Number)
  ChannelNumber <- Target |> pull(ChannelNumber)
  GainChannel <- Target |> pull(GainChannel)
  Max <- Target |> pull(Max)
  centerWaveLength <- TheWave
  bandWidth <- TheBand
  breakdown <- Target |> pull(breakdown)

  Output <- sprintf(
  '  <Detector Name="%s" Number="%s" ChannelNumber="%s" GainChannel="%s" Max="%s" centerWaveLength="%s" bandWidth="%s" breakdown="%s"/>'
  , Name, Number, ChannelNumber, GainChannel, Max, centerWaveLength, bandWidth, breakdown)
  
  return(Output)
}