#' Dashboard Internal, updates archive ApplicationLog.csv from Setup
#'
#' @param MainFolder The file.path to the Main Folder
#' @param x The Cytometer Folder Name
#'
#' @importFrom purrr map
#' @importFrom dplyr bind_rows arrange desc
#' @importFrom utils read.csv write.csv
#' @importFrom lubridate ymd_hms
#' @importFrom generics setdiff
#'
#' @return Updated tracking data CSV in the Archive Folder
#' @noRd
AppQCParse <- function(MainFolder, x) {

  Folder <- file.path(MainFolder, x)
  DailyQCFiles <- list.files(Folder, pattern="Application",
                                full.names = TRUE)

  if (!length(DailyQCFiles)==0) {

    if (length(DailyQCFiles)>=1) {
      Parsed <- map(.x=DailyQCFiles,
         .f=Luciernaga:::ApplicationLogParse) |> bind_rows()
      Parsed <- Parsed[grepl("cmd: SitFlush: Begin", Parsed$Command), ]
      Parsed <- Parsed |> arrange(desc(DateTime))
    } else {
      stop("Two csv files in the folder found!")
    }

    TheArchive <- file.path(Folder, "Archive")
    ArchivedDataFile <- list.files(TheArchive, pattern="Application",
                                   full.names = TRUE)

    if (!length(ArchivedDataFile)==0) {

      if (length(ArchivedDataFile)==1) {
        ArchivedData <- read.csv(ArchivedDataFile[1], check.names=FALSE)
      } else {
        message("Two csv files in the folder found!")
      }

      # Troubleshooting
      if (!ncol(ArchivedData) == ncol(Parsed)) {
        stop("Mismatched Number of Columns")
      }

      ArchivedData$DateTime <- ymd_hms(ArchivedData$DateTime)
      ArchivedData <- ArchivedData |> arrange(desc(DateTime))
      NewData <- generics::setdiff(Parsed, ArchivedData)
      UpdatedData <- rbind(NewData, ArchivedData)

      file.remove(ArchivedDataFile)

    } else {
      UpdatedData <- Parsed
    }

    file.remove(DailyQCFiles)

    UpdatedData <- UpdatedData |> arrange(desc(DateTime))

    name <- paste0("ApplicationData", x, ".csv")
    StorageLocation <- file.path(TheArchive, name)
    write.csv(UpdatedData, StorageLocation, row.names=FALSE)
  } else {
    message("No DailyQCFiles files to update with in ", x)
  }
}