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
  Temperature <- LaserSegment[TemperatureIndex]
  Temperature <- strsplit(Temperature, ",")[[1]]
  Temperature <- data.frame(Temperature=Temperature[[2]])
  Temperature$Temperature <- as.numeric(Temperature$Temperature)
  Temperature <- Temperature %>%
    mutate(across(everything(), ~ FALSE, .names = "Flag-{.col}"))

  LaserFrame <- LaserSegment[2:(FSCIndex-1)]
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
#' @param TrackChange Whether to derive a change between previous QC days.
#'
#' @importFrom purrr map
#' @importFrom dplyr select
#' @importFrom dplyr rename
#' @importFrom dplyr %>%
#' @importFrom tidyr starts_with
#' @importFrom lubridate mdy_hms
#' @importFrom lubridate mdy_hm
#'
#' @return A dataframe.
#' @export
#'
#' @examples
#' File_Location <- system.file("extdata", package = "Luciernaga")
#' CSV_Pattern <- ".CSV$"
#' CSV_Files <- list.files(path=File_Location, pattern=CSV_Pattern,
#'                        full.names=TRUE)
#' TidyData <- QC_FilePrep(CSV_Files, TrackChange = FALSE)

QC_FilePrep_LJTracking <- function(x, TrackChange, DailyQC){

    # Read the lines of the CSV file
    ReadInfo <- readLines(x)

    # More than one way into a house...
    ReadInfo <- ReadInfo[ReadInfo != ""]

    GainIndex <- grep("^Gain", ReadInfo)
    rCVIndex <- grep("^% rCV", ReadInfo)
    LaserIndex <- grep("^Laser Delay", ReadInfo)
    AreaIndex <- grep("^Area Scaling Factor", ReadInfo)
    FSCIndex <- grep("^FSC Area Scaling Factor", ReadInfo)
    LaserPowerIndex <- grep("^Laser Power", ReadInfo)

    ThePositions <- c(GainIndex, rCVIndex, LaserIndex, AreaIndex, FSCIndex,
                      LaserPowerIndex)
    #ThePositions <- c(GainIndex, LaserIndex, AreaIndex, rCVIndex,  FSCIndex,
    #LaserPowerIndex)

    ThePositions <- sort(ThePositions) #In case not everything was exported

    EmptyRowCheck <- logical(length(ThePositions))
    for (i in seq_along(ThePositions)) {
      if (ThePositions[i] > 1) {
        EmptyRowCheck[i] <- grepl("^,*$", ReadInfo[ThePositions[i] - 1])
      } else {
        EmptyRowCheck[i] <- FALSE
      }
    }

    #EmptyRowCheck

    StartPositions <- ThePositions + 1

    SpacePositions <- ThePositions - ifelse(EmptyRowCheck, 2, 0)
    NoSpacePositions <- SpacePositions - ifelse(EmptyRowCheck, 0, 1)

    EndPositions <- NoSpacePositions[-1]

    FinalPosition <- length(ReadInfo)

    #if (grepl("^,*$", ReadInfo[FinalPosition])){}

    ReadInfo[FinalPosition - 1]

    EndPositions <- c(EndPositions, FinalPosition)

    difference <- StartPositions - EndPositions

    if (length(unique(difference)) > 1) {
      stop("There are different number of blank rows in the .csv file.
           Please make sure   there is a single space in between chunks ")
    }

    TheLineChunks <- list()

    for (i in seq_along(StartPositions)) {
      start <- StartPositions[i]
      end <- EndPositions[i]

      TheLineChunks[[i]] <- ReadInfo[start:end]
    }

    for (i in seq_along(TheLineChunks)) {
      filename <- paste0("LineChunk_", i, ".txt")
      writeLines(TheLineChunks[[i]], filename)
      message("Lines extracted from ", StartPositions[i], " to ",
              EndPositions[i], " written to ", filename)
    }

    pattern <- "LineChunk_.*\\.txt$"
    #searchlocation <- dirname(path)
    TheLineChunkText <- list.files(pattern=pattern) #It's in working directory
    TheLineChunkText

    TidyedData <- map(TheLineChunkText, .f=ChunkReader) %>% bind_cols()

    unlink(TheLineChunkText) #For Cleanup (or make temp folder)

    TidyedData_Clean <- TidyedData[, colSums(is.na(TidyedData)) != nrow(TidyedData)]

    df <- TidyedData_Clean
    col_names <- colnames(df)

    for (i in seq_along(col_names)) {
      if (grepl("Out of Range Flag", col_names[i])) {
        prev_col_name <- col_names[i - 1]
        new_col_name <- paste0("Flag-", prev_col_name)
        names(df)[i] <- new_col_name
      }
    }

    RedundantColumns <- df %>% select(starts_with("DateTime")) %>%
      select(-1) %>% colnames()
    UpdatedDF <- df %>% select(-one_of(RedundantColumns)) %>%
      rename(DateTime = "DateTime...1")
    #colnames(UpdatedDF)
    #str(UpdatedDF)

    if (any(grepl("AM", UpdatedDF$DateTime))) {
      UpdatedDF$DateTime <- mdy_hms(UpdatedDF$DateTime)
    } else {UpdatedDF$DateTime <- mdy_hm(UpdatedDF$DateTime)}

    if (TrackChange == TRUE){
      FinalCol <- length(UpdatedDF)
      DateTime <- UpdatedDF %>% select(1)
      MainData <- UpdatedDF %>% select(all_of(2:FinalCol))

      Main_Columns <- MainData[, seq(1, ncol(MainData), by = 2)]
      Flag_Columns <- MainData[, seq(2, ncol(MainData), by = 2)]

      Main_Columns <- colnames(Main_Columns)
      Flag_Columns <- colnames(Flag_Columns)

      NewlyUpdatedDF <- map2(.x=Main_Columns, .y=Flag_Columns,
         .f=Internal_ChangeCalcs, TheData=MainData) %>% bind_cols()
      NewlyUpdatedDF <- cbind(DateTime, NewlyUpdatedDF)

    } else {NewlyUpdatedDF <- UpdatedDF}

    return(NewlyUpdatedDF)

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
