#' From Cytek Aurora Daily QC report, returns GainBaseline and RCVCutoffs
#'
#' @param x  Takes a Daily QC CSV file, and returns GainBaseline and RCVCutoffs template
#' @param outpath Default NULL, specifies location to store template .csv
#' @param returnType Options data or csv (saved to outpath)
#' 
#' @importFrom lubridate ymd
#' @importFrom lubridate hms
#' @importFrom dplyr mutate
#' @importFrom dplyr relocate
#' @importFrom dplyr select
#' @importFrom dplyr case_when
#' @importFrom dplyr rename
#' @importFrom stringr str_detect
#' @importFrom utils write.csv
#'
#' @return Either a data.frame or writes to a .csv file
#' 
#' @noRd
CytekDailyQC <- function(x, outpath=NULL, returnType){
  ReadInfo <- readLines(x)
  index <- grep("^Laser Settings", ReadInfo)
  Final <- length(ReadInfo)
  InitialIndex <- grep("^DailyQC", ReadInfo)
  Initial <- ReadInfo[InitialIndex]
  String <- gsub("DailyQCReport_", "", Initial)
  Parts <- strsplit(String, "_")[[1]]

  if (length(Parts) == 3){
      Instrument <- data.frame(Instrument = Parts[1])
      Date <- data.frame(Date = Parts[2])
      Date$Date <- ymd(Date$Date)
      Time <- data.frame(Time = Parts[3])
  } else if (length(Parts) == 2){
      Instrument <- "Unknown"
      Date <- data.frame(Date = Parts[1])
      Date$Date <- ymd(Date$Date)
      Time <- data.frame(Time = Parts[2])
  } else {stop("File Format for ", x, " not recognized")}

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

  if(any(colnames(TheData) == "DeltaGain")){
  TheData$DeltaGain <- as.integer(TheData$DeltaGain)
  }

  if(any(colnames(TheData) == "GainChange")){
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

  if (returnType == "csv"){
      if (is.null(outpath)){outpath <- getwd()}
      filename <- "CytekBaselineTemplate.csv"
      StorageLocation <- file.path(outpath, filename)
      write.csv(Updated, StorageLocation, row.names=FALSE)
  } else {return(Updated)}
}  

#' Imports Detector, GainBaseline and RCVCutoff for non-Cytek instruments
#'
#' @param x  A data.frame or file.path to the .csv with the Detector, GainBaseline and RCVCutoff
#'
#' @importFrom utils read.csv
#' 
#' @return Data for subsequent use
#' 
#' @noRd
NotCytekDailyQC <- function(x){
  if (!is.data.frame(x)){
      Data <- read.csv(x, check.names = FALSE)
  } else {Data <- x}

  TheseColumns <- c("Detector", "GainBaseline", "RCVCutoff")

  if (any(!TheseColumns %in% colnames(Data))){
      stop("CSV file should contain `Detector`, `GainBaseline` and `RCVCutoff` as column names")
  }

  Data$Detector <- as.character(Data$Detector)
  Data$GainBaseline <- as.numeric(Data$GainBaseline)
  Data$RCVCutoff <- as.numeric(Data$RCVCutoff)

  return(Data)
}


#' Converts the extracted BeadData into the ArchivedData format needed to plot Gain/RCV fails
#' 
#' @param data A data.frame object of the bead data
#' @param manufacturer Options Cytek or other, basically whether to use the DailyQC CSV or the template
#' @param baselinecutoffs The DailyQC .csv or the template .csv or associated file.path
#' @param returnTemplate Returns DailyQC template that can be adjusted beyond Cytek settings
#' @param outpath Default NULL, specifies where to returnTemplate .csv to
#' @param gainmultiplier Gain times this value is the cutoff point at which Gain Fails
#' 
#' @importFrom stringr str_detect
#' @importFrom dplyr mutate
#' @importFrom purrr map
#' @importFrom dplyr bind_cols
#' 
#' @return An updated data.frame containing the necessary Flag columns for plotting
#' 
#' @export
HolisticToArchived <- function(data, manufacturer="Cytek", baselinecutoffs, returnTemplate=FALSE, outpath=NULL, gainmultiplier=2){

    Internal <- colnames(data)[str_detect(colnames(data), "-A")]
    RCVs <- Internal[str_detect(Internal, "rCV")]
    Gains <- Internal[str_detect(Internal, "Gain")]
  
    if (manufacturer=="Cytek") {
      if (returnTemplate == TRUE){
        if (is.null(outpath)) {outpath <- getwd()}
        Cutoffs <- Luciernaga:::CytekDailyQC(x=baselinecutoffs,
        outpath=outpath, returnType="csv")
      } else {
        Cutoffs <- Luciernaga:::CytekDailyQC(x=baselinecutoffs, returnType="data")
        }
    } else {Cutoffs <- Luciernaga:::NotCytekDailyQC(x=baselinecutoffs)}
  
    Cutoffs <- Cutoffs |> mutate(GainBaseline=GainBaseline*gainmultiplier)
  
    TheRCVs <- map(.x=RCVs, .f=DerriveTheFlag, data=data, cutoffs=Cutoffs) |> bind_cols()
    TheGains <- map(.x=Gains, .f=DerriveTheFlag, data=data, cutoffs=Cutoffs) |> bind_cols()
    Assembled <- bind_cols(data, TheRCVs, TheGains)
    return(Assembled)
  }
  
#' Internal for Holistic to Archived, screens vs cutoff value, filling Flag column
#' 
#' @param x The column being iterated
#' @param data The data being selected on
#' @param cutoffs The processed Gain and RCV cutoff criteria
#' 
#' @importFrom stringr str_detect
#' @importFrom dplyr select
#' @importFrom dplyr filter
#' @importFrom dplyr case_when
#' @importFrom dplyr mutate
#' @importFrom dplyr pull
#' @importFrom tidyselect all_of
#' 
#' @return Individual flag columns
#' 
#' @noRd
DerriveTheFlag <- function(x, data, cutoffs){
  
    if (str_detect(x, "Gain")) {Type <- "Gain"
    } else if (str_detect(x, "rCV")) {Type <- "rCV"
    } else {Type <- "MFI"
    } 
  
    TheColumn <- gsub("-% rCV", "", gsub("_Gain", "", x))
    TheColumn <- gsub("-A", "", TheColumn)
    internaldata <- data |> dplyr::select(all_of(x))
    internalcutoffs <- Cutoffs |> dplyr::filter(Detector %in% TheColumn)
  
    if (nrow(internalcutoffs) == 1){
      if (Type == "Gain"){
        Comparison <- internalcutoffs |> dplyr::pull(GainBaseline)
      } else if (Type == "rCV") {
        Comparison <- internalcutoffs |> dplyr::pull(RCVCutoff)
      } else {message("Skipped")
      }
    }
  
    Status <- internaldata |> mutate(Flag = .data[[x]] >= Comparison)
    Status <- Status |> dplyr::select(-all_of(x))
    NewName <- paste0("Flag-", x)
    colnames(Status)[1] <- NewName
    return(Status)
  }
  