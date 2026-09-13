#' Dashboard Internal, loads updated data
#'
#' @param x The Cytometer Folder Name
#' @param MainFolder The file.path to the main folder
#' @param type Whether to return "MFI" or "Gain" plots
#'
#' @importFrom utils read.csv
#' @importFrom lubridate ymd_hms ymd hms mdy_hm
#' @importFrom stringr str_detect
#' @importFrom dplyr arrange desc
#'
#' @return Updated Tracking Data CSV for specified type
#' @noRd
CurrentData <- function(x, MainFolder, type) {

  ArchiveLocation <- file.path(MainFolder, x, "Archive")

  if (type == "MFI") {
    BeadData <- list.files(ArchiveLocation, pattern="Bead",
                           full.names=TRUE)
    Data <- read.csv(BeadData, check.names=FALSE)
    Data$DateTime <- ymd_hms(Data$DateTime)
    Data$DATE <- ymd(Data$DATE)
    Data$TIME <- hms(Data$TIME)
  }

  if (type == "Gain") {
    ArchiveData <- list.files(ArchiveLocation, pattern="Archived",
                              full.names=TRUE)
    Data <- read.csv(ArchiveData, check.names=FALSE)
    #lubridate::ymd_hms(Data$DateTime)

    if (any(str_detect(Data$DateTime, ":.*:"))) {
      Data$DateTime <- ymd_hms(Data$DateTime)
    } else {
      Data$DateTime <- mdy_hm(Data$DateTime)
    }
  }

  if (type == "Both") {
    BothData <- list.files(ArchiveLocation, pattern="Holistic",
                           full.names=TRUE)
    Data <- read.csv(BothData, check.names=FALSE)
    Data$DateTime <- ymd_hms(Data$DateTime)
    Data$DATE <- ymd(Data$DATE)
    Data$TIME <- hms(Data$TIME)
  }

  Data <- Data |> arrange(desc(DateTime))
  return(Data)
}