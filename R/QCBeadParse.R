#' Dashboard Internal, updates MFI tracking CSV
#'
#' @param x The cytometer folder name
#' @param MainFolder The file.path to main folder
#' @param timepointType Whether single or double.
#'
#' @importFrom flowWorkspace load_cytoset_from_fcs
#' @importFrom purrr map
#' @importFrom dplyr bind_rows mutate relocate arrange desc anti_join
#' @importFrom utils read.csv write.csv
#' @importFrom lubridate ymd_hms ymd hms
#'
#' @return Updated MFI tracking CSV
#' @noRd
QCBeadParse <- function(x, MainFolder, timepointType="double"){
  Folder <- file.path(MainFolder, x)
  FCS_Files <- list.files(Folder, pattern="fcs", full.names=TRUE)

  if(!length(FCS_Files) == 0){

    if (timepointType == "double"){
    QCBeads <- FCS_Files[grep("Before|After", FCS_Files)]
    } else {QCBeads <- FCS_Files}

    if (length(QCBeads) == 0){
      if (timepointType == "double"){
        message("No Before After detected in names")
      } else {message("No FCS Files detected")}
      QCBeads <- FCS_Files
    }

    BeforeAfter_CS <- tryCatch({
      load_cytoset_from_fcs(files = QCBeads, transformation = FALSE, truncate_max_range = FALSE)
    }, error = function(e) {
      Screen <- CytosetScreen(files = QCBeads)
      MainList <- which.max(sapply(Screen, length))
      Screen <- Screen[MainList][[1]]
      TheSave <- load_cytoset_from_fcs(files = Screen, transformation = FALSE, truncate_max_range = FALSE)

      message("Following error occurred: ", e$message, " Attempted Rescue by passing main list.",
      " Please run Luciernaga:::CytosetScreen to find skipped files in the smaller list entries")

      TheSave
    })

    if (timepointType == "double"){
      BeforeAfter <- map(.x=BeforeAfter_CS, .f=QC_GainMonitoring,
        sample.name = "TUBENAME", stats="median") |> bind_rows()
    } else {
      BeforeAfter <- map(.x=BeforeAfter_CS, .f=QC_GainMonitoring,
        sample.name = "$DATE", stats="median") |> bind_rows()
    }

    BeforeAfter <- BeforeAfter |> mutate(DateTime = DATE+TIME) |>
      relocate(DateTime, .before=DATE)

    BeforeAfter <- BeforeAfter |> arrange(desc(DateTime))

    ArchiveFolder <- file.path(Folder, "Archive")
    ArchiveCSV <- list.files(ArchiveFolder, pattern="Bead", full.names=TRUE)

    if (!length(ArchiveCSV) == 0){

    if (!length(ArchiveCSV) > 1){

      ArchiveData <- read.csv(ArchiveCSV, check.names=FALSE)
      ArchiveData$DateTime <- ymd_hms(ArchiveData$DateTime)
      ArchiveData$DATE <- ymd(ArchiveData$DATE)
      ArchiveData$TIME <- hms(ArchiveData$TIME)

      if (!ncol(BeforeAfter) == ncol(ArchiveData)){
        stop("Mismatched Number of Columns")
      }

      NewData <- BeforeAfter |>
        anti_join(ArchiveData, by = c("DATE", "TIME"))

      UpdatedData <- rbind(NewData, ArchiveData)

      file.remove(ArchiveCSV)

    } else {stop("Two BeadData csv files in the archive folder!")}

    } else {UpdatedData <- BeforeAfter}

    UpdatedData <- UpdatedData |> arrange(desc(DateTime))

    file.remove(FCS_Files)
    name <- paste0("BeadData", x, ".csv")
    StorageLocation <- file.path(ArchiveFolder, name)
    write.csv(UpdatedData, StorageLocation, row.names=FALSE)

    } else {message("No fcs files to update with in ", x)}
}