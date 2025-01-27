#' Parses Cytek Aurora Application Log for useful information
#' 
#' @param x A file.path to the individual ApplicationLog.txt file
#' @param returnType What to parse. Default is "clean" and returns the
#' operation data. Setting to "error" returns the Error and StackTrace notices
#' 
#' @importFrom lubridate mdy_hms
#' 
#' @return A data.frame with the selected returnType information
#' 
#' @noRd
ApplicationLogParse <- function(x, returnType="clean"){

  ReadInfo <- readLines(x)
  #ReadInfo <- ReadInfo[ReadInfo != ""]
  #Splits <- strsplit(ReadInfo, "\t")
  Splits <- strsplit(ReadInfo, "\\s*\\t\\s*")

  Data <- data.frame(
    DateTime = sapply(Splits, function(x) x[1]),
    Command = sapply(Splits, function(x) x[2])
  )
  
  Data$Command <- trimws(Data$Command, which="left")
  Data$Command <- gsub("----", "", Data$Command)
  
  Dataset <- ErrorStacks(data=Data, returnType=returnType)
  
  if (returnType != "clean"){
      Dataset <- Dataset[!is.na(Dataset$DateTime), ]
      return(Dataset)
  }
  
  Dataset <- Dataset[!is.na(Dataset$DateTime), ]

  if (nrow(Dataset) != 0){
  Dataset$DateTime <- trimws(Dataset$DateTime, which="right")
  Dataset$DateTime <- mdy_hms(Dataset$DateTime)
  #Troubleshooting <- Dataset[is.na(Dataset$DateTime), ]
  } else {message("No error-free rows detected, returning NULL, expect bind_rows to error")
    Dataset <- NULL}
  return(Dataset) 
  }

#' Internal for Application Log Parse, filters out or in the Error and StackTraces
#' 
#' @param data The passed data.frame
#' @param returnType The passed specification of data to keep, default is clean. 
#' 
#' @return The filtered clean or error data.frame
#' 
#' @noRd
ErrorStacks <- function(data, returnType="clean"){

    ErrorIndex <- which(grepl("^ERROR", data$Command))
    StackTraceIndex <- which(grepl("^STACKTRACE", data$Command))
    IndicatorChar <- substr(data$DateTime, 1, 1)
    ErrorReadoutsIndex <- which(grepl(" ", IndicatorChar))
    ErrorReadoutsIndex2 <- which(grepl("-", IndicatorChar))
    ErrorReadoutsIndex3 <- which(grepl("^[a-zA-Z]", IndicatorChar))
    ErrorReadoutsIndex4 <- which(grepl("'", IndicatorChar))
    
    These <- sort(unique(c(ErrorIndex, StackTraceIndex, ErrorReadoutsIndex,
     ErrorReadoutsIndex2, ErrorReadoutsIndex3, ErrorReadoutsIndex4)))
    
    if (returnType != "clean"){
        ErrorStacks <- data[These, ]
    } else {
        ErrorStacks <- data[-These, ]
        }
    return(ErrorStacks)
}
