#' Dashboard Internal, updates instrument tracking CSV
#'
#' @param MainFolder The file.path to the Main Folder
#' @param x The Cytometer Folder Name
#' @param Maintainer Logical override for when number of columns
#'   don't match
#'
#' @importFrom dplyr mutate across
#' @importFrom tidyselect starts_with
#' @importFrom utils read.csv write.csv
#' @importFrom lubridate ymd_hms
#' @importFrom generics setdiff
#'
#' @return Updated tracking data CSV in the Archive Folder
#' @noRd
LevyJenningsParse <- function(MainFolder, x, Maintainer=FALSE) {

  Folder <- file.path(MainFolder, x)
  LJTrackingFiles <- list.files(Folder, pattern="LJ",
                                full.names = TRUE)

  if (!length(LJTrackingFiles)==0) {

    if (length(LJTrackingFiles)==1) {
      Parsed <- QC_FilePrep(x=LJTrackingFiles, TrackChange=FALSE)
      Parsed <- Parsed |>
        mutate(across(starts_with("Flag"), ~ as.logical(.)))
    } else {
      stop("Two csv files in the folder found!")
    }

    TheArchive <- file.path(Folder, "Archive")
    ArchivedDataFile <- list.files(TheArchive, pattern="Archived",
                                   full.names = TRUE)

    if (!length(ArchivedDataFile)==0) {

      if (length(ArchivedDataFile)==1) {
        ArchivedData <- read.csv(ArchivedDataFile, check.names=FALSE)
      } else {
        message("Two csv files in the folder found!")
      }

      # Troubleshooting
      if (!ncol(ArchivedData) == ncol(Parsed)) {

        if (Maintainer==TRUE) {
          TheseColumns <- setdiff(colnames(Parsed), colnames(ArchivedData))

          for (col in TheseColumns) {
            ArchivedData[[col]] <- NA
          }

          if (!ncol(ArchivedData) == ncol(Parsed)) {
            stop("Still no rescue")
          }
        } else {
          stop("The number of columns for the new data don't match
       that of the archived data. Please make sure to
       export the Levy-Jennings trackings with all available
       parameters")
        }
      }

      ArchivedData$DateTime <- ymd_hms(ArchivedData$DateTime)
      ArchivedData <- ArchivedData |>
        mutate(across(starts_with("Flag"), ~ as.logical(.)))
      NewData <- generics::setdiff(Parsed, ArchivedData)
      UpdatedData <- rbind(NewData, ArchivedData)

      file.remove(ArchivedDataFile)

    } else {
      UpdatedData <- Parsed
    }

    file.remove(LJTrackingFiles)

    name <- paste0("ArchivedData", x, ".csv")
    StorageLocation <- file.path(TheArchive, name)
    write.csv(UpdatedData, StorageLocation, row.names=FALSE)
  } else {
    message("No LevyJennings files to update with in ", x)
  }
}