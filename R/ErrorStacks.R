#' Internal for Application Log Parse, filters out or in the
#'  Error and StackTraces
#' 
#' @param data The passed data.frame
#' @param returnType The passed specification of data to keep,
#'  default is clean. 
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