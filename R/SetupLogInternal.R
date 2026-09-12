
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