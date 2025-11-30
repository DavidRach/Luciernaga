#' Forgive me Cytek, for I hath sinned, this shit and giggles
#' function attempts to generate a config file to use to unmix other
#' manufacturers files with SpectroFlo. It also serves as a learning 
#' tool for what elements in the config file do what. Original nerd-snipe
#' credit goes to the Cytometry Discord discussion
#' 
#' @param NumberDetectors QC_ReferenceLibrary input to specify alternate instrument
#' 
#' @importFrom dplyr pull bind_rows left_join mutate case_when
#' @importFrom purrr map
#' @importFrom stringr str_detect
#' 
#' @noRd
FrankensteinsConfig <- function(NumberDetectors){

Data <- Luciernaga:::InstrumentReferences(NumberDetectors) |> pull(Detector) |> unique()
AllDetectors <- length(Data)
  
if (any(str_detect(Data, "-"))){
  Lasers <- sub("-.*", "", Data) |> unique()
} else {Lasers <- str_remove_all(Data, "[0-9]") |> unique()}
  
# x <- Lasers[2] 
# data <- Data
DetectorCount <- map(.x=Lasers, .f=DetectorsPerLaser, data=Data) |> bind_rows()

TargetLasers <- data.frame(Lasers=Lasers, check.names=FALSE)
#str(TargetLasers)

TargetLasers <- left_join(TargetLasers, DetectorCount, by="Lasers")

ReferenceChart <- TargetLasers |> mutate(Original = case_when(
  Lasers == "320" ~ "355",
  Lasers == "349" ~ "355",
  Lasers == "355" ~ "355",
  Lasers == "405" ~ "405",
  Lasers == "445" ~ "405",
  Lasers == "488" ~ "488",
  Lasers == "561" ~ "561",
  Lasers == "637" ~ "640",
  Lasers == "640" ~ "640",
  Lasers == "785" ~ "640",
  Lasers == "808" ~ "640",
  Lasers == "UV" ~ "355",
  Lasers == "V" ~ "405",
  Lasers == "B" ~ "488",
  Lasers == "YG" ~ "561",
  Lasers == "R" ~ "640",
  TRUE ~ "")
)

TheseUnits <- ReferenceChart |> pull(Original) |> unique()

Summary <- ReferenceChart |> mutate(ID="") |> mutate(Name="") |>
  mutate(Wavelength="") |> mutate(Position="") |> 
  mutate(Delay="") |> mutate(IsReference="") |> 
  mutate(GainControlBoard="") |>
  mutate(Start="") |> mutate(End="")

Summary <- Summary |> mutate(ID = case_when(
  Original == "355" ~ "104",
  Original == "405" ~ "100",
  Original == "488" ~ "101",
  Original == "561" ~ "103",
  Original == "640" ~ "102",
  TRUE ~ "")
)

Summary <- Summary |> mutate(Name = case_when(
  Original == "355" ~ "UV",
  Original == "405" ~ "Violet",
  Original == "488" ~ "Blue",
  Original == "561" ~ "YellowGreen",
  Original == "640" ~ "Red",
  TRUE ~ "")
)

Summary <- Summary |> mutate(Wavelength = case_when(
  Original == "355" ~ "355",
  Original == "405" ~ "405",
  Original == "488" ~ "488",
  Original == "561" ~ "561",
  Original == "640" ~ "640",
  TRUE ~ "")
)

Summary <- Summary |> mutate(Position = case_when(
  Original == "355" ~ "4",
  Original == "405" ~ "1",
  Original == "488" ~ "2",
  Original == "561" ~ "0",
  Original == "640" ~ "3",
  TRUE ~ "")
)

Summary <- Summary |> mutate(Delay = case_when(
  Original == "355" ~ "40",
  Original == "405" ~ "-20",
  Original == "488" ~ "0",
  Original == "561" ~ "-40",
  Original == "640" ~ "20",
  TRUE ~ "")
)

Summary <- Summary |> mutate(GainControlBoard = case_when(
  Original == "355" ~ "6-7",
  Original == "405" ~ "0-1",
  Original == "488" ~ "2-3",
  Original == "561" ~ "2-5",
  Original == "640" ~ "4-5",
  TRUE ~ "")
)

Summary <- Summary |> mutate(IsReference = case_when(
  Original == "355" ~ "False",
  Original == "405" ~ "False",
  Original == "488" ~ "True",
  Original == "561" ~ "False",
  Original == "640" ~ "False",
  TRUE ~ "")
)
  
Metadata <- FranksHeadGears(NumberDetectors = NumberDetectors)
#writeLines(Metadata, "Metadata.xml")
  
CorrectOrder <- c("561", "405", "488", "640", "355")
TheseUnits <- CorrectOrder[CorrectOrder %in% TheseUnits]

# x <- TheseUnits[3]
Combined <- map(.f=FranksArms, .x=TheseUnits, reference=Summary,
  data=Data)
  
AllLasers <- paste(Combined, collapse = "\n  ")
FinalOutput <- gsub("#PlaceLasersHere", AllLasers, Metadata, fixed = TRUE)
writeLines(FinalOutput, "Final.xml")
return(FinalOutput)
}

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

#' Internal FrankensteinsConfig,returns overall .xml config
#' 
#' @param NumberDetectors Number detectors of alternate instrument
#' 
#' @noRd
FranksHeadGears <- function(NumberDetectors){
  InstrumentAndDetectors <- paste("Frankenstein", NumberDetectors, sep=" ")

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


#' Internal FrankensteinsConfig, determines number of detectors per laser
#' 
#' @param x The iterated laser
#' @param data The vector of detectors
#' 
#' @importFrom stringr str_detect
#' 
#' @noRd
DetectorsPerLaser <- function(x, data){
  Lasers <- x
  if (Lasers == "V"){
    Lasers <- paste0("^", Lasers)
    Conditional <- TRUE
  } else {Conditional <- FALSE}

  DetectorCount <- length(data[str_detect(data, Lasers)])

  if (Conditional == TRUE){Lasers <- "V"}

  Stats <- data.frame(Lasers, DetectorCount, check.names=FALSE)
  return(Stats)
}
