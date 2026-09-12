#' Converts Cytek Aurora Daily QC report into a data frame.
#'
#' @param x  Takes a Daily QC CSV file, and converts into a "tidyed" dataframe for plotting.
#' Currently works on our 3L, 4L, 5L Auroras.
#'
#' @importFrom lubridate ymd hms
#' @importFrom dplyr mutate relocate select case_when rename across
#' @importFrom tidyr pivot_wider
#' @importFrom tidyselect everything
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