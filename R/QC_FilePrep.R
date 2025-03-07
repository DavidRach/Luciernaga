#' Converts Cytek Aurora Daily QC report into a data frame.
#'
#' @param x  Takes a Daily QC CSV file, and converts into a "tidyed" dataframe for plotting.
#' Currently works on our 3L, 4L, 5L Auroras.
#'
#' @importFrom lubridate ymd
#' @importFrom lubridate hms
#' @importFrom dplyr mutate
#' @importFrom dplyr relocate
#' @importFrom dplyr select
#' @importFrom dplyr case_when
#' @importFrom dplyr rename
#' @importFrom tidyr pivot_wider
#' @importFrom tidyselect everything
#' @importFrom dplyr across
#'
#' @return A dataframe.
#' @export
#'
#' @examples
#' File_Location <- system.file("extdata", package = "Luciernaga")
#' CSV_Pattern <- ".CSV$"
#' CSV_Files <- list.files(path=File_Location, pattern=CSV_Pattern,
#'                        full.names=TRUE)
#' TidyData <- QC_FilePrep_DailyQC(CSV_Files)

QC_FilePrep_DailyQC <- function(x){
  ReadInfo <- readLines(x)
  #ReadInfo
  index <- grep("^Laser Settings", ReadInfo)
  Final <- length(ReadInfo)

  #Instrument
  InitialIndex <- grep("^DailyQC", ReadInfo)
  Initial <- ReadInfo[InitialIndex]
  String <- gsub("DailyQCReport_", "", Initial)
  Parts <- strsplit(String, "_")[[1]]

  if (length(Parts) == 3){
    Instrument <- data.frame(Instrument = Parts[1])
    Date <- data.frame(Date = Parts[2])
    Date$Date <- lubridate::ymd(Date$Date)
    Time <- data.frame(Time = Parts[3])
  } else if (length(Parts) == 2){
      Instrument <- "Unknown"
      Date <- data.frame(Date = Parts[1])
      Date$Date <- lubridate::ymd(Date$Date)
      Time <- data.frame(Time = Parts[2])
  } else {stop("File Format for ", x, " not recognized")}
  
  Time$Time <- sub("(\\d{2})(\\d{2})(\\d{2})", "\\1:\\2:\\3", Time$Time)
  Time$Time <- lubridate::hms(Time$Time)
  Intro <- cbind(Date, Time, Instrument)
  Intro <- Intro |>
    mutate(DateTime=Date+Time) |>
    relocate(DateTime, .before=1) |>
    select(-Date, -Time)

  # Detector Section
  StartDetector <- grep("^Laser,Detector", ReadInfo)
  #FinalDetector <- grep("^Red,R8", ReadInfo)
  DetectorSegment <- ReadInfo[StartDetector:(index-3)]
  DetectorLength <- length(DetectorSegment)

  header <- strsplit(DetectorSegment[1], ",")[[1]]
  #extra <- strsplit(DetectorSegment[2], ",")[[1]]
  #header <- c(header, extra)

  data <- DetectorSegment[3:DetectorLength]
  data <- strsplit(data, ",")

  TheData <- do.call(rbind, lapply(
    data, function(x) as.data.frame(t(x), stringsAsFactors = FALSE)))
  colnames(TheData) <- header

  TheData$Detector <- gsub(" .*", "", TheData$Detector)
  colnames(TheData) <- gsub(" ", "", colnames(TheData))
  TheData$Gain <- as.integer(TheData$Gain)
  TheData$`%rCV` <- as.numeric(TheData$`%rCV`)

  if(any(colnames(TheData) == "DeltaGain")){
  TheData$DeltaGain <- as.integer(TheData$DeltaGain)
  }

  if(any(colnames(TheData) == "GainChange")){
    TheData <- TheData %>% rename(DeltaGain = GainChange)
    TheData$DeltaGain <- as.integer(TheData$DeltaGain)
  }


  #TheData[2,7] <- 7

  Updated <- TheData %>%
    mutate(Baseline=Gain-DeltaGain) %>%
    mutate(Comparison=Baseline*2) %>%
    mutate(GainFlag = Gain > Comparison) %>%
    mutate(RCVFlag = `%rCV` > 6) %>%
    mutate(RCVFlag = case_when(
      Detector == "SSC" | Detector == "SSC-B" & `%rCV` < 8 ~ FALSE,
      TRUE ~ RCVFlag))

  TheGain <- Updated %>% select(Detector, Gain)
  TheGain$Detector <- paste0(TheGain$Detector, "-Gain")
  TheGain <- TheGain %>%
    pivot_wider(., names_from = Detector, values_from = Gain)

  GainFlagged <- Updated %>% select(Detector, GainFlag)
  GainFlagged$Detector <- paste0("Flag-", GainFlagged$Detector, "-Gain")
  GainFlagged <- GainFlagged %>%
    pivot_wider(., names_from = Detector, values_from = GainFlag)

  TheRCV <- Updated %>% select(Detector, `%rCV`)
  TheRCV$Detector <- paste0(TheRCV$Detector, "-% rCV")
  TheRCV <- TheRCV %>%
    pivot_wider(., names_from = Detector, values_from = `%rCV`)

  RCVFlagged <- Updated %>% select(Detector, RCVFlag)
  RCVFlagged$Detector <- paste0("Flag-", RCVFlagged$Detector, "-% rCV")
  RCVFlagged <- RCVFlagged %>%
    pivot_wider(., names_from = Detector, values_from = RCVFlag)

  Assembly <- cbind(Intro, TheGain, GainFlagged, TheRCV, RCVFlagged)

  # Laser Section
  LaserSegment <- ReadInfo[index:Final]

  FSCIndex <- grep("^FSC", LaserSegment)
  FSCAreaScaling <- LaserSegment[FSCIndex]
  FSCAreaScaling <- strsplit(FSCAreaScaling, ",")[[1]]
  FSCArea = data.frame(FSCAreaScalingFactor=FSCAreaScaling[[2]])
  FSCArea$FSCAreaScalingFactor <- as.numeric(FSCArea$FSCAreaScalingFactor)
  FSCArea <- FSCArea %>%
    mutate(across(everything(), ~ FALSE, .names = "Flag-{.col}"))


  FlowRateIndex <- grep("^Flow", LaserSegment)

  if(length(FlowRateIndex) > 0){
  FlowRate <- LaserSegment[FlowRateIndex]
  FlowRate <- strsplit(FlowRate, ",")[[1]]
  FlowRate = data.frame(FlowRate=FlowRate[[2]])
  FlowRate$FlowRate <- as.numeric(FlowRate$FlowRate)
  FlowRate <- FlowRate %>%
    mutate(across(everything(), ~ FALSE, .names = "Flag-{.col}"))
  FSCArea <- cbind(FSCArea, FlowRate)
  }

  TemperatureIndex <- grep("^Temperature", LaserSegment)

  if (length(TemperatureIndex) > 0){
    Temperature <- LaserSegment[TemperatureIndex]
    Temperature <- strsplit(Temperature, ",")[[1]]
    Temperature <- data.frame(Temperature=Temperature[[2]])
    Temperature$Temperature <- as.numeric(Temperature$Temperature)
    Temperature <- Temperature %>%
      mutate(across(everything(), ~ FALSE, .names = "Flag-{.col}"))
  } else {
    Temperature <- data.frame(Temperature=0, `Flag-Temperature`=FALSE, check.names = FALSE)
  }

  if (length(FSCIndex) == 2){
  LaserFrame <- LaserSegment[2:(FSCIndex[1]-1)]
  } else {LaserFrame <- LaserSegment[2:(FSCIndex-1)]}
    
  LaserFrame <- strsplit(LaserFrame, ",")

  LaserFrame <- do.call(rbind, lapply(
    LaserFrame, function(x) as.data.frame(t(x), stringsAsFactors = FALSE)))

  colnames(LaserFrame) <- LaserFrame[1,]
  LaserFrame <- LaserFrame[-1,]

  colnames(LaserFrame) <- gsub(" L", "L", colnames(LaserFrame))
  colnames(LaserFrame) <- gsub(" A", "A", colnames(LaserFrame))
  #LaserFrame <- LaserFrame %>% mutate(across(everything(), as.numeric))

  if (any(str_detect(colnames(LaserFrame), "Laser Power"))){

    LaserPower <- LaserFrame %>%
      select(-`Area Scaling Factor`, -`Laser Delay`) %>%
      pivot_wider(
        names_from = Laser,
        values_from = `Laser Power`,
        names_glue = "{Laser}-Laser Power"
      )

    LaserPower <- LaserPower %>% mutate(across(everything(), as.numeric))
    LaserPower <- LaserPower %>%
      mutate(across(everything(), ~ FALSE, .names = "Flag-{.col}"))

    LaserFrame <- LaserFrame %>% select(-`Laser Power`)

    FSCArea <- cbind(LaserPower, FSCArea)
  }

  LaserDelay <- LaserFrame %>%
    select(-`Area Scaling Factor`) %>%
    pivot_wider(
      names_from = Laser,
      values_from = `Laser Delay`,
      names_glue = "{Laser}-Laser Delay"
    )

  LaserDelay <- LaserDelay %>% mutate(across(everything(), as.numeric))
  LaserDelay <- LaserDelay %>%
    mutate(across(everything(), ~ FALSE, .names = "Flag-{.col}"))

  LaserArea <- LaserFrame %>%
    select(-`Laser Delay`) %>%
    pivot_wider(
      names_from = Laser,
      values_from = `Area Scaling Factor`,
      names_glue = "{Laser}-Area Scaling Factor"
    )

  LaserArea <- LaserArea %>% mutate(across(everything(), as.numeric))
  LaserArea <- LaserArea %>%
    mutate(across(everything(), ~ FALSE, .names = "Flag-{.col}"))

  NewData <- cbind(Assembly, LaserDelay, LaserArea, FSCArea, Temperature)

  NewData <- NewData %>% select(-Instrument)

  return(NewData)
  }

#' Converts the Cytek Aurora (TM)'s QC report into a data frame.
#'
#' @param x  Takes a Levy-Jennings QC tracking report saved as a .csv file, and
#' converts into a "tidyed" dataframe for plotting. Currently works on our 3L, 4L, 5L
#' Auroras. Please reach out if you find an issue, the .csv export varies a bit and
#' I want to continue to improve on the code to handle these odd exceptions.
#' @param DailyQC A single DailyQCReport .csv file, used to import baseline settings.
#'
#' @importFrom purrr map2
#' @importFrom tidyr pivot_longer
#' @importFrom tidyr pivot_wider
#' @importFrom tidyselect all_of
#' @importFrom dplyr "%>%"
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @importFrom dplyr left_join
#' @importFrom dplyr case_when
#' @importFrom dplyr bind_cols
#' @importFrom dplyr across
#' @importFrom tidyselect everything
#' @importFrom lubridate mdy_hms
#'
#' @return A dataframe.
#' @export
#'
#' @examples
#' File_Location <- system.file("extdata", package = "Luciernaga")
#' CSV_Pattern <- ".CSV$"
#' CSV_Files <- list.files(path=File_Location, pattern=CSV_Pattern,
#'                        full.names=TRUE)
#' TidyData <- QC_FilePrep(x=CSV_Files, DailyQC=DailyQC)

QC_FilePrep_LJTracking <- function(x, DailyQC){
  ReadInfo <- readLines(x)
  ReadInfo <- ReadInfo[ReadInfo != ""]

  GainIndex <- grep("^Gain", ReadInfo)
  rCVIndex <- grep("^% rCV", ReadInfo)
  LaserIndex <- grep("^Laser Delay", ReadInfo)
  AreaIndex <- grep("^Area Scaling Factor", ReadInfo)
  FSCIndex <- grep("^FSC Area Scaling Factor", ReadInfo)
  LaserPowerIndex <- grep("^Laser Power", ReadInfo)

  ThePositions <- c(GainIndex, rCVIndex, LaserIndex,
                    AreaIndex, FSCIndex, LaserPowerIndex)

  Final <- length(ReadInfo)
  StartPositions <- ThePositions + 1
  EndPositions <- ThePositions - 1
  EndPositions <- EndPositions[-1]
  EndPosition <- c(EndPositions, Final)

  if (StartPositions[length(StartPositions)] > Final){
    StartHere <- StartPositions[-length(StartPositions)]
    EndHere <- EndPosition[-length(EndPosition)]
  } else {
    StartHere <- StartPositions
    EndHere <- EndPosition
  }

  TheFrames <- map2(.x=StartHere, .y=EndHere, ReadInfo=ReadInfo, .f=ParseThis)

  Updated <- TheFrames

  Gains <- Updated[[1]]
  GainLength <- ncol(Gains)

  GainLonger <- Gains %>%
    pivot_longer(all_of(2:GainLength), names_to="DetectorMain", values_to="Value")

  Baseline <- DailyQCBaseline(DailyQC)
  Comparison <- cbind(GainLonger, Baseline)
  Comparison$Value <- as.numeric(Comparison$Value)
  Comparison$Comparison <- as.numeric(Comparison$Comparison)

  Comparison <- Comparison %>%
    mutate(GainFlag = Value > Comparison) %>%
    select(-Detector, -Comparison) %>% mutate(DetectorFlag=paste0("Flag-", DetectorMain))

  GainsMain <- Comparison %>% select(DateTime, DetectorMain, Value) %>%
    pivot_wider(names_from = DetectorMain, values_from = Value)

  GainsFlag <- Comparison %>% select(DateTime, DetectorFlag, GainFlag) %>%
    pivot_wider(names_from = DetectorFlag, values_from = GainFlag)

  TheGains <- left_join(GainsMain, GainsFlag, by="DateTime")

  RCVs <- Updated[[2]]
  RCVsLength <- ncol(RCVs)
  RCVsLonger <- RCVs %>%
    pivot_longer(all_of(2:RCVsLength), names_to = "RCVMain", values_to = "Value")

  RCVsLonger$Value <- as.numeric(RCVsLonger$Value)

  Comparison <- RCVsLonger %>% mutate(RCVFlag = Value > 6) %>%
    mutate(Detector=gsub("-% rCV", "", RCVMain)) %>%
    mutate(RCVFlag = case_when(
      Detector == "SSC" | Detector == "SSC-B" & Value < 8 ~ FALSE,
      TRUE ~ RCVFlag)) %>% select(-Detector) %>%
      mutate(FlagRCV=paste0("Flag-", RCVMain))

  RCVsMain <- Comparison %>% select(DateTime, RCVMain, Value) %>%
    pivot_wider(names_from=RCVMain, values_from = Value)

  RCVsFlag <- Comparison %>% select(DateTime, FlagRCV, RCVFlag) %>%
    pivot_wider(names_from=FlagRCV, values_from = RCVFlag)

  TheRCVs <- left_join(RCVsMain, RCVsFlag, by="DateTime")

  UpdatedItems <- length(Updated)

  Others <- Updated[3:UpdatedItems] %>%
    lapply(function(df) df %>% select(-DateTime)) %>%
    bind_cols()

  Others[] <- lapply(Others, as.numeric)

  OtherAdditional <- Others %>%
    mutate(across(everything(), ~ FALSE, .names = "Flag-{.col}"))

  TheDataset <- left_join(TheGains, TheRCVs, by="DateTime")
  TheDataset <- cbind(TheDataset, OtherAdditional)

  TheDataset$DateTime <- mdy_hms(TheDataset$DateTime)

  TheDataset$DateTime <-
  return(TheDataset)
}

#' Internal for LJTracking parser
#'
#' @param x A DailyQC report .csv
#'
#' @importFrom lubridate ymd
#' @importFrom lubridate hms
#' @importFrom dplyr mutate
#' @importFrom dplyr relocate
#' @importFrom dplyr select
#'
#' @return A data.frame comparing the cutoffs for Gains
#' @noRd
DailyQCBaseline <- function(x){
  ReadInfo <- readLines(x)
  #ReadInfo
  index <- grep("^Laser Settings", ReadInfo)
  Final <- length(ReadInfo)

  #Instrument
  InitialIndex <- grep("^DailyQC", ReadInfo)
  Initial <- ReadInfo[InitialIndex]
  String <- gsub("DailyQCReport_", "", Initial)
  Parts <- strsplit(String, "_")[[1]]
  Instrument <- data.frame(Instrument = Parts[1])
  Date <- data.frame(Date = Parts[2])
  Date$Date <- lubridate::ymd(Date$Date)
  Time <- data.frame(Time = Parts[3])
  Time$Time <- sub("(\\d{2})(\\d{2})(\\d{2})", "\\1:\\2:\\3", Time$Time)
  Time$Time <- lubridate::hms(Time$Time)
  Intro <- cbind(Date, Time, Instrument)
  Intro <- Intro %>%
    mutate(DateTime=Date+Time) %>%
    relocate(DateTime, .before=1) %>%
    select(-Date, -Time)

  # Detector Section
  StartDetector <- grep("^Laser,Detector", ReadInfo)
  #FinalDetector <- grep("^Red,R8", ReadInfo)
  DetectorSegment <- ReadInfo[StartDetector:(index-3)]
  DetectorLength <- length(DetectorSegment)

  header <- strsplit(DetectorSegment[1], ",")[[1]]
  #extra <- strsplit(DetectorSegment[2], ",")[[1]]
  #header <- c(header, extra)

  data <- DetectorSegment[3:DetectorLength]
  data <- strsplit(data, ",")

  TheData <- do.call(rbind, lapply(
    data, function(x) as.data.frame(t(x), stringsAsFactors = FALSE)))
  colnames(TheData) <- header

  TheData$Detector <- gsub(" .*", "", TheData$Detector)
  colnames(TheData) <- gsub(" ", "", colnames(TheData))
  TheData$Gain <- as.integer(TheData$Gain)
  TheData$DeltaGain <- as.integer(TheData$DeltaGain)
  TheData$`%rCV` <- as.numeric(TheData$`%rCV`)

  #TheData[2,7] <- 7

  Updated <- TheData %>%
    mutate(Baseline=Gain-DeltaGain) %>%
    mutate(Comparison=Baseline*2)

  Cutoff <- Updated %>% select(Detector, Comparison)
  return(Cutoff)
}

#' Internal for LJTracking Parse
#'
#' @param x The starting index
#' @param y The End Index
#' @param ReadInfo The ReadLines object
#'
#' @importFrom stringr str_detect
#' @importFrom dplyr select
#' @importFrom tidyselect all_of
#'
#' @return A data.frame chunk from the original
#'
#' @noRd
ParseThis <- function(x, y, ReadInfo){
  DetectorSegment <- ReadInfo[x:y]
  DetectorLength <- length(DetectorSegment)
  header <- strsplit(DetectorSegment[1], ",")[[1]]
  data <- strsplit(DetectorSegment, ",")
  TheData <- do.call(rbind, lapply(data, function(x) as.data.frame(t(x), stringsAsFactors = FALSE)))
  colnames(TheData) <- header
  TheData <- TheData[-1,]
  cols_to_select <- which(str_detect(colnames(TheData), "Out of Range Flag"))
  SelectedData <- TheData %>% select(-all_of(cols_to_select))
  return(SelectedData)
}

#' Internal for CytekQCFilePrep
#'
#' @importFrom utils read.csv
#'
#' @return An internal value
#'
#' @noRd
ChunkReader <- function(x){
  ReadChunks <- read.csv(x, check.names = FALSE)
}


#' Internal for CytekQCFilePrep
#'
#' @importFrom purrr map
#' @importFrom dplyr select
#' @importFrom dplyr rename
#' @importFrom tidyr starts_with
#' @importFrom lubridate mdy_hms
#' @importFrom lubridate mdy_hm
#'
#' @return An internal value
#'
#' @noRd
Internal_ChangeCalcs <- function(x, y, TheData){
  xx <- TheData %>% select(all_of(x))
  yy <- TheData %>% select(all_of(y))
  FinalValue <- xx[nrow(xx),]
  z <- xx[2:nrow(xx),]
  z <- c(z, FinalValue)
  Prelim <- xx %>% mutate(TheSubtract = z)
  Result <- Prelim %>% mutate(Crazy = .[[1]] - .[[2]])
  FinalResult <- Result %>% select(-2)
  FinalResult <- FinalResult %>% bind_cols(yy)
  FinalResult <- FinalResult %>% mutate(FlagCrazy = .[[3]])
  DiffColName <- paste0("Change_", colnames(FinalResult)[1])
  DiffFlagColName <- gsub("Flag-", "Flag-Change_", colnames(FinalResult)[3])
  colnames(FinalResult)[2] <- DiffColName
  colnames(FinalResult)[4] <- DiffFlagColName
  #FinalResult

  return(FinalResult)
}
