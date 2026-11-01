#' Dashboard Internal, loads updated data
#'
#' @param x The Cytometer Folder Name
#' @param MainFolder The file.path to the main folder
#' @param type Whether to return "MFI" or "Gain" plots
#'
#' @importFrom utils read.csv
#' @importFrom lubridate ymd_hms ymd mdy hms mdy_hm parse_date_time
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
    if (grepl("^\\d{4}[-/]\\d{1,2}[-/]\\d{1,2}$", x)) {Data$DATE <- ymd(Data$DATE)
    } else if (grepl("^\\d{1,2}[-/]\\d{1,2}[-/]\\d{4}$", x)) {Data$DATE <- mdy(Data$DATE)}
    Data$TIME <- hms(Data$TIME)
  }

  if (all(is.na(Data$DateTime))){Data <- Data |>
    mutate(DateTime = parse_date_time(paste(DATE, TIME),
                                      orders = c("mdy HMS", "ymd HMS")))}

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