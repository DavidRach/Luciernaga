#' Takes the folder outputs from, parses the stain index for all the
#' combinations, returning as a large long data.frame
#' 
#' @param folder_location File.path to the parent folder where all the subfolders containing
#' the variant .fcs files were stored
#' @param outpath Where you want to store the processed data
#' @param excludeThese Used to exclude columns from transformation, default is "FSC|SSC|Time"
#' @param channelRange Default for biexponential transformation is 4096
#' @param maxValue Default for biexponential transformation is 4194304
#' @param pos Default for biexponential transformation is 5.62
#' @param neg Default for biexponential transformation is 0
#' @param widthBasis Default for biexponential transformation is -1000
#' @param inverse.transform Default is TRUE
#' @param stringAppend Default is -A
#' 
#' @importFrom purrr map
#' @importFrom dplyr bind_rows
#' @importFrom utils write.csv
#' 
#' @return A data.frame containing staining index for all the FMO folders fcs files. 
#' 
#' @export
#' 
#' @examples A <- 2 + 2
#' 
StainBrightnessIndexCalculator <- function(folder_location,
                                            outpath = NULL,
                                            excludeThese = "FSC|SSC|Time",
                                            channelRange = 4096,
                                            maxValue = 4194304,
                                            pos = 5.62,
                                            neg = 0,
                                            widthBasis = -1000,
                                            inverse.transform = TRUE,
                                            stringAppend = "-A") {

  AllUnmix <- list.files(folder_location, full.names = TRUE,
                          pattern = "AllUnmix")
  FMO <- list.files(folder_location, full.names = TRUE, pattern = "FMO")
  FMO_Folders <- list.files(FMO, full.names = TRUE)
  SingleUnmix <- list.files(folder_location, full.names = TRUE,
                             pattern = "SingleUnmix")

  # Handle the FMOs
  # x <- FMO_Folders[26]
  TheData <- purrr::map(.x = FMO_Folders,
                         BrightnessIndexIterator,
                         excludeThese = excludeThese,
                         channelRange = channelRange,
                         maxValue = maxValue,
                         pos = pos,
                         neg = neg,
                         widthBasis = widthBasis,
                         inverse.transform = inverse.transform,
                         stringAppend = stringAppend,
                         .progress = TRUE)

  TheData <- TheData |> bind_rows()

  if (is.null(outpath)) {
    outpath <- getwd()
  }

  filename <- "FMO_StainIndex.csv"
  StoreHere <- file.path(outpath, filename)
  write.csv(TheData, StoreHere, row.names = FALSE)

  # Handle the 'Full-Unmixed'

  return(TheData)
}