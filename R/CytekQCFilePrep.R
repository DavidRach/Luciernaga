#' Converts the Cytek Aurora (TM)'s QC report into a data frame.
#'
#' @param x  Takes a Levy-Jennings QC tracking report saved as a .csv file, and converts into a "tidyed" dataframe for plotting.
#'           Currently works on our 3L, 4L, 5L Auroras. Please reach out if you find an issue, the .csv export varies a bit and
#'           I want to continue to improve on the code to handle these odd exceptions.
#' @param TrackChange Whether to derive a change between previous QC days.
#'
#' @importFrom purrr map
#' @importFrom dplyr select
#' @importFrom dplyr rename
#' @importFrom tidyr starts_with
#' @importFrom lubridate mdy_hms
#' @importFrom lubridate mdy_hm
#'
#' @return A dataframe.
#' @export
#'
#' @examples NULL

CytekQCFilePrep <- function(x, TrackChange){

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

    ThePositions <- c(GainIndex, rCVIndex, LaserIndex, AreaIndex, FSCIndex, LaserPowerIndex)
    #ThePositions <- c(GainIndex, LaserIndex, AreaIndex, rCVIndex,  FSCIndex, LaserPowerIndex)

    ThePositions <- sort(ThePositions) #In case not everything was exported

    EmptyRowCheck <- logical(length(ThePositions))
    for (i in seq_along(ThePositions)) {
      if (ThePositions[i] > 1) {
        EmptyRowCheck[i] <- grepl("^,*$", ReadInfo[ThePositions[i] - 1])
      } else {
        EmptyRowCheck[i] <- FALSE
      }
    }

    EmptyRowCheck

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
      stop("There are different number of blank rows in the .csv file. Please make sure   there is a single space in between chunks ")
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
      cat("Lines extracted from", StartPositions[i], "to", EndPositions[i], "written to", filename, "\n")
    }

    pattern <- "LineChunk_.*\\.txt$"
    #searchlocation <- dirname(path)
    TheLineChunkText <- list.files(pattern=pattern) #It's in working directory
    TheLineChunkText

    .ChunkReader <- function(x){
      ReadChunks <- read.csv(x, check.names = FALSE)
    }

    TidyedData <- map(TheLineChunkText, .f=.ChunkReader) %>% bind_cols()

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

    RedundantColumns <- df %>% select(starts_with("DateTime")) %>% select(-1) %>% colnames()
    UpdatedDF <- df %>% select(-one_of(RedundantColumns)) %>% rename(DateTime = "DateTime...1")
    #colnames(UpdatedDF)
    #str(UpdatedDF)

    if (any(grepl("AM", UpdatedDF$DateTime))) {UpdatedDF$DateTime <- mdy_hms(UpdatedDF$DateTime)
    } else {UpdatedDF$DateTime <- mdy_hm(UpdatedDF$DateTime)}

    if (TrackChange == TRUE){
      FinalCol <- length(UpdatedDF)
      DateTime <- UpdatedDF %>% select(1)
      MainData <- UpdatedDF %>% select(all_of(2:FinalCol))

      Main_Columns <- MainData[, seq(1, ncol(MainData), by = 2)]
      Flag_Columns <- MainData[, seq(2, ncol(MainData), by = 2)]

      Main_Columns <- colnames(Main_Columns)
      Flag_Columns <- colnames(Flag_Columns)

      .Internal_ChangeCalcs <- function(x, y, TheData){
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

      NewlyUpdatedDF <- map2(.x=Main_Columns, .y=Flag_Columns, .f=.Internal_ChangeCalcs, TheData=MainData) %>% bind_cols()
      NewlyUpdatedDF <- cbind(DateTime, NewlyUpdatedDF)

    } else {NewlyUpdatedDF <- UpdatedDF}

    return(NewlyUpdatedDF)

  }
