#' Internal for LJTracking parser
#'
#' @param x A DailyQC report .csv
#'
#' @importFrom lubridate ymd hms
#' @importFrom dplyr mutate relocate select
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