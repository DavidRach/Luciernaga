#' Internal FrankensteinsConfig, works on each Laser .xml
#' 
#' @param x The laser being iterated on
#' @param data The vector of detectors
#' @param reference the Summary metadata
#' 
#' @importFrom dplyr filter pull mutate row_number
#' @importFrom stringr str_detect str_replace fixed
#' @importFrom purrr map
#' 
#' @noRd
FranksArms <- function(x, data, reference){

TheTarget <- reference |> filter(Original %in% x)
TheseLasers <- TheTarget |> pull(Lasers)
  
if (nrow(TheTarget) > 2){
  Detectors <- data[str_detect(data, paste(TheseLasers, collapse = "|"))]
} else {
  Detectors <- data[str_detect(data, paste(TheseLasers, collapse = "|"))]
}
  
TotalDetectors <- length(Detectors) 

ID <- TheTarget |> pull(ID)
Name <- TheTarget |> pull(Name)
Wavelength <- TheTarget |> pull(Wavelength)
Position <- TheTarget |> pull(Position)
Delay <- TheTarget |> pull(Delay)
IsReference <- TheTarget |> pull(IsReference)

Boards <- TheTarget |> pull(GainControlBoard)
Boards <- strsplit(Boards, "-")[[1]]
Board1 <- Boards[1]
Board2 <- Boards[2]

if(!TheTarget$Original == "561"){
    BoardChunk <- sprintf(
    '<Board Number="%s" HeaderSize="8" DataSize="64">
    </Board>
    <Board Number="%s" HeaderSize="12" DataSize="80">
    #PlaceDetectorsHere
    </Board>', Board2, Board1)
    TheBoard <- Board1
} else {
    BoardChunk <- sprintf(
    '<Board Number="%s" HeaderSize="8" DataSize="64">
    #PlaceDetectorsHere  
    </Board>
    <Board Number="%s" HeaderSize="12" DataSize="80">
    </Board>', Board2, Board1)
   TheBoard <- Board2
}

  LaserChunk <- sprintf(
  '<Laser
      Id="%s"
      Name="%s"
      WaveLength="%s"
      Position="%s"
      Delay="%s"
      GainControlBoard="%s"   
      IsReference="%s">
    #PlaceBoardsHere
  </Laser>', ID, Name, Wavelength, Position, Delay, TheBoard, IsReference)

  Combined <- str_replace(LaserChunk, fixed("#PlaceBoardsHere"), BoardChunk)
  #writeLines(Combined, "Lasers.xml")
  #return(Combined)

  #Now to derrive the detector chunk. 
  Detectors <- Detectors[length(Detectors):1]
  Detectors <- data.frame(Detectors=Detectors, check.names=FALSE)
  Detectors <- Detectors |> mutate(Name="") |> mutate(Number=row_number()) |> mutate(ChannelNumber=row_number()) |>
    mutate(GainChannel=row_number()-1) |> mutate(Max="10000") |> mutate(centerWavelength="") |>
    mutate(bandWidth="") |> mutate(breakdown="150")

  Iterators <- Detectors |> pull(Detectors)

  # x <- Iterators[1]
  # detector <- Detectors
  DetectorOutputs <- map(.x=Iterators, detector=Detectors, .f=FranksToes)
  AllDetectors <- paste(DetectorOutputs, collapse = "\n    ")
  LaserOutput <- gsub("#PlaceDetectorsHere", AllDetectors, Combined, fixed = TRUE)
  return(LaserOutput)
}