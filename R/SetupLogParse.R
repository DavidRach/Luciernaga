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
#' @importFrom lubridate mdy_hms
#' @importFrom dplyr filter
#' @importFrom stringr str_starts
#' @importFrom purrr map
#' @importFrom dplyr bind_rows
#' @importFrom dplyr mutate
#' @importFrom lubridate floor_date
#' @importFrom lubridate ymd_hms
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



#' Internal for SetupEngineLogLibraryParse
#' 
#' @param data The dataframe from which ranges to filter
#' @param x A list containing ranges to filter.  
#' 
#' @importFrom dplyr select
#' @importFrom stringr str_extract
#' @importFrom tibble tibble
#' 
#' @return A dataframe row
#' 
#' @noRd
SetupLogInternal <- function(x, data){
  
  if (!is.list(x)){x <- list(x)
    }

  InternalData <- data[x[[1]],]
  Date <- InternalData[1,1]
  Name <- InternalData[1,2]
  Name <- gsub("Ref Control Name ", "", Name)
  Reference <- InternalData[3,2]

  Subset <- InternalData[-1:-3,] |> select(Comment)

  if (nrow(Subset) < 2) {
    stop("Not enough rows in Subset to process. Skipping this range.")
    return(NULL)
  }

  Main <- Subset[seq(1, nrow(Subset), by = 2), ]
  SOV <- Subset[seq(2, nrow(Subset), by = 2), ]

  Detector <- str_extract(Main, "^[^:]+")
  Positive_MFI <- str_extract(Main,
   "(?<=Positive MFI = )\\d+\\.\\d+") |> as.numeric()
  Negative_MFI <- str_extract(Main,
   "(?<=Negative MFI = )\\d+\\.\\d+") |> as.numeric()
  Normalized <- str_extract(SOV,
   "(?<=Sov = )[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?") |> as.numeric()
  
  Dataset <- tibble(Detector=Detector, PositiveMFI=Positive_MFI,
   NegativeMFI=Negative_MFI, Normalized=Normalized)

  Dataset <- data.frame(Date, Name, Reference, Dataset)
  return(Dataset)
}


#' A wrapper for taking Reference Library Signatures derrived from SetupLog and generate
#' plots for each of them
#' 
#' @param x A
#' @param distinguish A
#' @param data A
#' @param columnname A
#' @param detectorcolumn A
#' @param valuecolumn A
#' @param Normalize A
#' @param TheFormat A
#' 
#' @importFrom dplyr filter
#' @importFrom dplyr pull
#' @importFrom purrr map
#' @importFrom purrr flatten
#' 
#' @return A ggplot2 object
#' 
#' @noRd 
SetupLibraryPlotWrapper <- function(x, distinguish, data, columnname, detectorcolumn, 
  valuecolumn, Normalize, TheFormat){
  
  SubsetData <- data |> filter(.data[[columnname]] %in% x)
  Distinguisher <- SubsetData |> pull(distinguish) |> unique()

  ThePlots <- map(.x=Distinguisher, .f=InternalSetupLibraryWrap, distinguish=distinguish,
   data=SubsetData, columnname = columnname, detectorcolumn=detectorcolumn,
    valuecolumn=valuecolumn, Normalize=Normalize, TheFormat=TheFormat)

  ThePlots <- flatten(ThePlots)

  return(ThePlots)
}

#' Internal for SetupLibraryPlotWrapper
#' 
#' @param x A
#' @param distinguish A
#' @param data A
#' @param columnname A
#' @param detectorcolumn A
#' @param valuecolumn A
#' @param Normalize A
#' @param TheFormat A
#' 
#' @importFrom dplyr filter
#' @importFrom dplyr pull
#' @importFrom purrr map
#' 
#' @return A ggplot2 object
#' 
#' @noRd
InternalSetupLibraryWrap <- function(x, distinguish, data, columnname,
 detectorcolumn, valuecolumn, Normalize, TheFormat){

  SubsetData <- data |> filter(.data[[distinguish]] %in% x)
  TheTarget <- SubsetData |> pull(columnname) |> unique()

  InnerPlots <- map(.x=TheTarget, .f=QC_ViewSignature, data=SubsetData,
   columnname = columnname, detectorcolumn=detectorcolumn,
    valuecolumn=valuecolumn, Normalize=Normalize, TheFormat=TheFormat)

   return(InnerPlots)
}