#' From Cytek Aurora Daily QC report, returns GainBaseline and RCVCutoffs
#'
#' @param x Takes a Daily QC CSV file, and returns GainBaseline and
#'   RCVCutoffs template
#' @param outpath Default NULL, specifies location to store template .csv
#' @param returnType Options data or csv (saved to outpath)
#'
#' @importFrom lubridate ymd hms
#' @importFrom dplyr mutate relocate select case_when rename
#' @importFrom stringr str_detect
#' @importFrom utils write.csv
#'
#' @return Either a data.frame or writes to a .csv file
#'
#' @noRd
CytekDailyQC <- function(x, outpath=NULL, returnType="data") {

  if (length(x) > 1) {
    ReadInfo <- x
  } else {
    ReadInfo <- readLines(x)
  }
  index <- grep("^Laser Settings", ReadInfo)
  Final <- length(ReadInfo)
  InitialIndex <- grep("^DailyQC", ReadInfo)
  Initial <- ReadInfo[InitialIndex]
  String <- gsub("DailyQCReport_", "", Initial)
  Parts <- strsplit(String, "_")[[1]]

  if (length(Parts) == 3) {
    Instrument <- data.frame(Instrument = Parts[1])
    Date <- data.frame(Date = Parts[2])
    Date$Date <- ymd(Date$Date)
    Time <- data.frame(Time = Parts[3])
  } else if (length(Parts) == 2) {
    Instrument <- "Unknown"
    Date <- data.frame(Date = Parts[1])
    Date$Date <- ymd(Date$Date)
    Time <- data.frame(Time = Parts[2])
  } else {
    stop("File Format for ", x, " not recognized")
  }

  Time$Time <- sub("(\\d{2})(\\d{2})(\\d{2})", "\\1:\\2:\\3", Time$Time)
  Time$Time <- hms(Time$Time)
  Intro <- cbind(Date, Time, Instrument)
  Intro <- Intro |>
    mutate(DateTime=Date+Time) |>
    relocate(DateTime, .before=1) |>
    select(-Date, -Time)

  StartDetector <- grep("^Laser,Detector", ReadInfo)
  DetectorSegment <- ReadInfo[StartDetector:(index-3)]
  DetectorLength <- length(DetectorSegment)
  header <- strsplit(DetectorSegment[1], ",")[[1]]
  data <- DetectorSegment[3:DetectorLength]
  data <- strsplit(data, ",")

  TheData <- do.call(rbind, lapply(
      data, function(x) as.data.frame(t(x), stringsAsFactors = FALSE)))
  colnames(TheData) <- header

  TheData$Detector <- gsub(" .*", "", TheData$Detector)
  colnames(TheData) <- gsub(" ", "", colnames(TheData))
  TheData$Gain <- as.integer(TheData$Gain)
  TheData$`%rCV` <- as.numeric(TheData$`%rCV`)

  if (any(colnames(TheData) == "DeltaGain")) {
    TheData$DeltaGain <- as.integer(TheData$DeltaGain)
  }

  if (any(colnames(TheData) == "GainChange")) {
    TheData <- TheData |> rename(DeltaGain = GainChange)
    TheData$DeltaGain <- as.integer(TheData$DeltaGain)
  }

  Updated <- TheData |>
    mutate(GainBaseline=Gain-DeltaGain) |>
    select(Detector, GainBaseline) |> mutate(RCVCutoff=NA)
  Updated$RCVCutoff <- as.numeric(Updated$RCVCutoff)
  Updated <- Updated |> mutate(RCVCutoff = case_when(
    str_detect(Detector, "SSC") ~ 8,
    TRUE ~ 6
  ))

  Updated$Detector <- as.character(Updated$Detector)
  Updated$GainBaseline <- as.numeric(Updated$GainBaseline)
  Updated$RCVCutoff <- as.numeric(Updated$RCVCutoff)

  if (returnType == "csv") {
    if (is.null(outpath)) {
      outpath <- getwd()
    }
    filename <- "CytekBaselineTemplate.csv"
    StorageLocation <- file.path(outpath, filename)
    write.csv(Updated, StorageLocation, row.names=FALSE)
  } else {
    return(Updated)
  }
}