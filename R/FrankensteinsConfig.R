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