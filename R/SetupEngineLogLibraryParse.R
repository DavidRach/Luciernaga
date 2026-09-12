#' Processes SetupEngineLog.csv to a tidy format, and returns Reference Library signatures found within. 
#' 
#' @param x The file.path to the desired SetupEngineLog.csv
#' @param returnArg Default is "Ref", can be adjusted any value when returnType = "list" to parse other
#' portions of the tidyed log
#' @param returnType Default is data to retrieve the Reference Library signatures, alternative is list to
#' return the tidyed data without filtering or processing
#' @param NumberDetectors The Aurora number of detectors, used to retrieve the complete information
#'  rather cutting off abruptly
#'
#' @importFrom lubridate mdy_hms floor_date ymd_hms
#' @importFrom dplyr filter bind_rows mutate
#' @importFrom stringr str_starts
#' @importFrom purrr map
#' 
#' @return A data.frame object
#' 
#' @noRd 
SetupEngineLogLibraryParse <- function(x, returnArg="Ref",
 returnType="data", NumberDetectors=64){

  ReadInfo <- readLines(x)
  HashLines <- grep("^#", ReadInfo)
  Starts <- HashLines[seq(1, length(HashLines), by = 2)]
  Stops <- HashLines[seq(2, length(HashLines), by = 2)]
  ToRemove <- unlist(Map(function(start, end) seq(start, end), Starts, Stops))
  NoHash <- ReadInfo[-ToRemove]
  NotEmpty <- NoHash[NoHash != ""]

  Split <- strsplit(NotEmpty, "\t")
  Data <- do.call(rbind, lapply(Split, function(x) {
    data.frame(DateTime = x[1], Comment = x[2], stringsAsFactors = FALSE)}))
  Data$DateTime <- mdy_hms(Data$DateTime)
  MissedASpot <- is.na(Data$DateTime)
  Data <- Data[!MissedASpot, ]

  if(returnType == "List"){
  Hmm <- Data |> filter(str_starts(Comment, returnArg))
  return(Hmm)
  } else {
    Hmm <- which(str_starts(Data$Comment, returnArg))
    UpTo <- (NumberDetectors*2)+3

    Ranges <- lapply(Hmm, function(start){
      end <- start + UpTo - 1
      if (end > nrow(Data)) {return(NULL)}
      seq(start, end)
    })
    #Handles any overshooting at end
    Ranges <- Filter(Negate(is.null), Ranges) 

    # Conventional Modes don't contain "Measurement A values"
    RealRanges <- Ranges[sapply(Ranges, function(Verify) {RowCheck <- Verify[2]
      if (RowCheck <= nrow(Data)) {return(Data$Comment[RowCheck] == "Measurement A")
        } else {return(FALSE)}})]
    
    Intermediate <- map(.x=RealRanges, .f=SetupLogInternal, data=Data) |> bind_rows()
    
    # Attempt to facilitate duplicate removal subsequently
    Final <- Intermediate |> mutate(Date = floor_date(ymd_hms(Date), unit = "5 minutes"))
    return(Final)
  }
}