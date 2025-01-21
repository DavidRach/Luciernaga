#' Creates a Description List for .fcs file from scratch.
#'
#' @param x The data.frame of exprs values (including time and scatters)
#'
#' @importFrom purrr map
#' @importFrom purrr flatten
#' @importFrom stringr str_detect
#' @importFrom dplyr filter
#' @importFrom dplyr pull
#' @importFrom dplyr mutate
#' @importFrom dplyr relocate
#' @importFrom dplyr select
#' @importFrom dplyr case_when
#'
#' @return A TBD product
#' @noRd
DescriptionGenesis <- function(x){
 Parameters <- ParameterPrep(x)

 Description <- list(
   FCSversion = '3',
   '$BEGINANALYSIS' = '0',
   '$BEGINDATA' = '45270',
   '$BEGINSTEXT' = '0',
   '$BTIM' = '04:00:00.00',
   '$BYTEORD' = '4,3,2,1',
   '$CYT' = 'Aurora',
   '$CYTOLIB_VERSION' = '2.16.0',
   '$CYTSN' = 'U010524',
   '$DATATYPE' = 'F',
   '$DATE' = '05-Jan-2024',
   '$ENDANALYSIS' = '0',
   '$ENDDATA' = '3005269',
   '$ENDSTEXT' = '0',
   '$ETIM' = '04:01:00.00',
   '$FIL' = 'Example.fcs',
   '$FLOWRATE' = 'Medium',
   '$INST' = '0',
   '$MODE' = 'L',
   '$NEXTDATA' = '0',
   '$OP' = 'FulanoDeTal'
 )
 #View(Description)

 ParameterAdjust <- Parameters
 ParameterAdjust <- RowNameReorder(ParameterAdjust)
 NewParameters <- rownames(ParameterAdjust)
 TheParamList <- map(.x=NewParameters, .f=ParameterDictate, data=ParameterAdjust)
 TheParamList <- flatten(TheParamList)
 #View(TheParamList)

 TheNames <- Parameters[['name']]
 TheNames <- TheNames[!grepl("Time|FSC|SSC", TheNames, ignore.case = TRUE)]
 Diagonally <- diag(1, nrow = length(TheNames), ncol = length(TheNames))
 colnames(Diagonally) <- TheNames
 Diagonally[2, 1] <- 0.000001

 # SpillOver Matrix Creation

 Description_2 <- list(
   '$PAR'='74',
   '$PROJ'='2025_Example',
   '$SPILLOVER'= Diagonally,
   '$TIMESTEP'='0.0001',
   '$TOT'='10000',
   '$VOL'='30.28',
   'APPLY COMPENSATION'='FALSE',
   'CHARSET'='utf-8',
   'CREATOR'='Luciernaga 0.99.1',
   'FILENAME'='C:\\Users\\FulanoDeTal\\Desktop\\ExampleData.fcs',
   'FSC ASF'='1.04',
   'GROUPNAME'='ExampleGroup',
   'GUID'='ExampleData.fcs'
 )

 #View(Description_2)

 Lasers <- list()
 if(any(str_detect(Parameters$name, "UV"))){Lasers <- c(Lasers, "UV")
 }
 if(any(str_detect(Parameters$name, "^V"))){Lasers <- c(Lasers, "V")
 }
 if(any(str_detect(Parameters$name, "^B"))){Lasers <- c(Lasers, "B")
 }
 if(any(str_detect(Parameters$name, "^YG"))){Lasers <- c(Lasers, "YG")
 }
 if(any(str_detect(Parameters$name, "^R"))){Lasers <- c(Lasers, "R")
 }
 Lasers <- unlist(Lasers)

 LaserData <- data.frame(LaserOrder=c("YG", "V", "B", "R", "UV"),
                         ASF=c('1.11', '1.1', '1.03', '1.13', '1.19'),
                         DELAY=c('-40.875', '-20.65', '0', '19.2', '38.975'),
                         NAME=c("YellowGreen", "Violet",
                                     "Blue", "Red", "UV"))
 LaserData_subset <- LaserData %>% dplyr::filter(LaserOrder %in% Lasers) %>%
   mutate(LaserNumber=row_number()) %>% relocate(LaserNumber, .before=ASF)
 Lasers <- LaserData_subset %>% pull(LaserOrder)
 AllLasers <- map(.x=Lasers, .f=LaserList, data=LaserData_subset)
 AllLasers <- flatten(AllLasers)
 #View(AllLasers)

 DisplaySetup <- ParameterAdjust %>% select(name)
 DisplaySetup <- DisplaySetup %>% mutate(Display="LOG")
 DisplaySetup <- DisplaySetup %>% mutate(
   Display = case_when(name == "Time" ~ 'LIN', TRUE ~ Display),
   Display = case_when(name == "FSC-A" ~ 'LIN', TRUE ~ Display),
   Display = case_when(name == "FSC-H" ~ 'LIN', TRUE ~ Display),
   Display = case_when(name == "FSC-W" ~ 'LIN', TRUE ~ Display),
   Display = case_when(name == "SSC-A" ~ 'LIN', TRUE ~ Display),
   Display = case_when(name == "SSC-H" ~ 'LIN', TRUE ~ Display),
   Display = case_when(name == "SSC-W" ~ 'LIN', TRUE ~ Display),
   Display = case_when(name == "SSC-B-A" ~ 'LIN', TRUE ~ Display),
   Display = case_when(name == "SSC-B-H" ~ 'LIN', TRUE ~ Display),
   Display = case_when(name == "SSC-B-W" ~ 'LIN', TRUE ~ Display))

 TheDisplayNames <- rownames(DisplaySetup)

 DisplayList <- map(.x=TheDisplayNames, .f=DisplayInternal, data=DisplaySetup)
 Display <- flatten(DisplayList)

 #View(Display)

 Description_3 <- list(
   'THRESHOLD'='(FSC,150000)',
   'TUBENAME'='ExampleData',
   'USERSETTINGNAME'='2025_Simulation',
   'WINDOW EXTENSION'='3',
   'ORIGINALGUID'='ExampleData.fcs'
 )
 #View(Description_3)

 ###################
 # List, Assemble! #
 ###################

 TheMegaList <- c(Description, TheParamList, Description_2,
                  AllLasers, Display, Description_3)

 return(TheMegaList)
}


#' Internal for Chorizo, setup Display parameters for Description
#'
#' @param x The iterated rowname
#' @param data The Display data.frame to be filtered
#'
#' @importFrom dplyr pull
#'
#' @return The list of display parameters for the detector
#' @noRd
DisplayInternal <- function(x, data){
  Subset <- data[x,]

  TheX <- gsub("$", "", fixed=TRUE, x)

  Display <- paste0(TheX, "DISPLAY")
  DisplayVal <- Subset %>% pull(Display)

  DisplayList <- list(Display=DisplayVal)
  names(DisplayList) <- Display
  return(DisplayList)
}

#' Internal for Chorizo, iterates out individual detector parameters
#'
#' @param x Iterated In Dollar-P-Number for filtering
#' @param data The reordered parameter data
#'
#' @importFrom dplyr pull
#' @importFrom stringr str_detect
#'
#' @return The individual detectors list of parameters for description file
#' @noRd
ParameterDictate <- function(x, data){
  Subset <- data[x,]

  Bits <- paste0(x, "B")
  Ehh <- paste0(x, "E")
  Name <-  paste0(x, "N")
  Range <-  paste0(x, "R")
  Type <-  paste0(x, "TYPE")
  Volt <-  paste0(x, "V")

  BitsVal <- '32'
  EhhVal <- '0,0'
  NameVal <- Subset %>% pull(name)
  RangeVal <- Subset %>% pull(range)
  RangeVal <- as.character(RangeVal)

  if(str_detect(NameVal, "FSC")){
    TypeVal <-  'Forward_Scatter'
  } else if (str_detect(NameVal, "SSC")){
    TypeVal <-  'Side_Scatter'
  } else if (str_detect(NameVal, "Time")){
    TypeVal <-  'Time'
  } else {TypeVal <-  'Raw_Fluorescence'}

  if (!str_detect(NameVal, "Time")){
  VoltVal <-  '307'
  DetectorList <- list(
    Bits = BitsVal,
    Ehh = EhhVal,
    Name = NameVal,
    Range = RangeVal,
    Type = TypeVal,
    Volt = VoltVal
  )
  names(DetectorList) <- c(Bits, Ehh, Name, Range, Type, Volt)
  } else {
    DetectorList <- list(
      Bits = BitsVal,
      Ehh = EhhVal,
      Name = NameVal,
      Range = RangeVal,
      Type = TypeVal
    )
    names(DetectorList) <- c(Bits, Ehh, Name, Range, Type)
  }

  return(DetectorList)
}

#' Internal for Chorizo, rearranges dollar-p-number according to bizarre order
#'
#' @param x The parameter data.frame
#'
#' @return The rearranged parameter data.frame according to wacky order.
#' @noRd
RowNameReorder <- function(x){
  TheRowNames <- rownames(x)
  LastElement <- TheRowNames[length(TheRowNames)]
  LastElement <- gsub("$P", "", fixed=TRUE, LastElement)
  LastNumber <- as.numeric(LastElement)

  TheIntegers <- 1:LastNumber

  TheOrder <- SplitThemUp(TheIntegers, RangeStart=10, RangeSize=10)

  Rearranged <- x[TheOrder, ]

  return(Rearranged)
}

#' Internal for Chorizo, determines the Parameter row name order
#'
#' @param TheIntegers The iterated number of row names
#' @param RangeStart Default 10, because it matches
#' @param RangeSize Default 10, because it matches
#'
#' @importFrom purrr flatten
#'
#' @return A vector of dollar-p-number names to rearrange parameter data
#'
#' @noRd
SplitThemUp <- function(TheIntegers, RangeStart, RangeSize) {
  TheIntegers <- TheIntegers[TheIntegers >= RangeStart]
  RangesList <- split(TheIntegers, (TheIntegers - RangeStart) %/% RangeSize + 1)
  RangesLength <- length(RangesList)
  InitialAppend <- 1:RangesLength
  Remainder <- 1:9
  These <- setdiff(Remainder, InitialAppend)

  for (i in seq_along(InitialAppend)) {
    RangesList[[i]] <- c(RangesList[[i]], InitialAppend[i])
  }

  RangesList[[RangesLength]] <- c(RangesList[[RangesLength]], These)

  RangesList <- flatten(RangesList)
  TheList <- unlist(RangesList)
  TheList <- paste0("$P", TheList)

  return(TheList)
}

#' Internal for Chorizo, returns Laser parameters for Description
#'
#' @param x The iterated laser
#' @param data The filtered laser list
#'
#' @importFrom dplyr filter
#' @importFrom dplyr pull
#'
#' @return The list of laser parameters for the Description list
#'
#' @noRd
LaserList <- function(x, data){
  Subset <- data %>% dplyr::filter(LaserOrder %in% x)
  Number <- Subset %>% pull(LaserNumber)

  LaserASF <- paste0("LASER", Number, "ASF")
  LaserDelay <- paste0("LASER", Number, "DELAY")
  LaserName <- paste0("LASER", Number, "NAME")

  LaserASFVal <- Subset %>% pull(ASF)
  LaserDelayVal <- Subset %>% pull(DELAY)
  LaserNameVal <- Subset %>% pull(NAME)

  LaserList <- list(
    LaserASF = LaserASFVal,
    LaserDelay = LaserDelayVal,
    LaserName = LaserNameVal
  )

  # Optionally, set the names programmatically if dynamic
  names(LaserList) <- c(LaserASF, LaserDelay, LaserName)

  return(LaserList)
}

#' Internal for Chorizo, generates a parameter data.frame from provided data.
#'
#' @param x A data.frame with parameter column names in the correct order.
#'
#' @importFrom dplyr mutate
#' @importFrom dplyr pull
#' @importFrom dplyr case_when
#' @importFrom dplyr row_number
#'
#' @return A data.frame with the parameter data
#' @noRd
ParameterPrep <- function(x){
  TheNames <- colnames(x)
  ParameterStandin <- data.frame(name=TheNames, desc=NA, range=4194303,
                                 minRange=-111.00000, maxRange=4194303)
  TheParams <- ParameterStandin %>%
    mutate(TheParams=paste0("$P", row_number())) %>% pull(TheParams)
  rownames(ParameterStandin) <- TheParams
  ParameterStandin$desc <- as.character(ParameterStandin$desc)

  ParameterStandin <- ParameterStandin %>% mutate(
    range = case_when(name == "Time" ~ 532116, TRUE ~ range),
    minRange = case_when(name == "Time" ~ 0.00000	, TRUE ~ minRange),
    minRange = case_when(name == "SSC-W" ~ 0.00000	, TRUE ~ minRange),
    minRange = case_when(name == "SSC-H" ~ 0.00000	, TRUE ~ minRange),
    minRange = case_when(name == "SSC-A" ~ 0.00000	, TRUE ~ minRange),
    minRange = case_when(name == "FSC-W" ~ 0.00000	, TRUE ~ minRange),
    minRange = case_when(name == "FSC-H" ~ 0.00000	, TRUE ~ minRange),
    minRange = case_when(name == "FSC-A" ~ 0.00000	, TRUE ~ minRange),
    minRange = case_when(name == "SSC-B-W" ~ 0.00000	, TRUE ~ minRange),
    minRange = case_when(name == "SSC-B-H" ~ 0.00000	, TRUE ~ minRange),
    minRange = case_when(name == "SSC-B-A" ~ 0.00000	, TRUE ~ minRange),
    minRange = case_when(name == "UV8-A" ~ -91.25813	, TRUE ~ minRange),
    minRange = case_when(name == "V6-A" ~ 0.00000	, TRUE ~ minRange),
    minRange = case_when(name == "V7-A" ~ 0.00000	, TRUE ~ minRange),
    minRange = case_when(name == "V8-A" ~ 0.00000	, TRUE ~ minRange),
    minRange = case_when(name == "V9-A" ~ 0.00000	, TRUE ~ minRange),
    minRange = case_when(name == "V10-A" ~ 0.00000	, TRUE ~ minRange),
    minRange = case_when(name == "B3-A" ~ 0.00000	, TRUE ~ minRange)
  )

  return(ParameterStandin)
}

