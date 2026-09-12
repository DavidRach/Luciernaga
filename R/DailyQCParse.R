#' Dashboard Internal, updates instrument tracking CSV from DailyQC
#'
#' @param MainFolder The file.path to the Main Folder
#' @param x The Cytometer Folder Name
#'
#' @importFrom purrr map
#' @importFrom dplyr bind_rows
#' @importFrom dplyr mutate
#' @importFrom dplyr across
#' @importFrom tidyselect starts_with
#' @importFrom utils read.csv
#' @importFrom lubridate ymd_hms
#' @importFrom generics setdiff
#' @importFrom dplyr arrange
#' @importFrom dplyr desc
#' @importFrom utils write.csv
#' @importFrom lubridate ymd
#'
#' @return Updated tracking data CSV in the Archive Folder
#' @noRd
DailyQCParse <- function(MainFolder, x){

  Folder <- file.path(MainFolder, x)
  DailyQCFiles <- list.files(Folder, pattern="DailyQC",
                                full.names = TRUE)

  if (!length(DailyQCFiles)==0){

    if (length(DailyQCFiles)>=1){

      Parsed <- map(.x=DailyQCFiles, .f=QC_FilePrep_DailyQC) |> bind_rows()
      Parsed <- Parsed |> mutate(across(starts_with("Flag"), ~ as.logical(.)))

    } else {stop("Two csv files in the folder found!")}
    
    # New Integration # Verify that it adds correctly
      ShinyData <- ShinyQCSummary(x=Parsed, Instrument=x)
      HistoricalPath <- file.path(MainFolder, "HistoricalData.csv")
      History <- list.files(MainFolder, pattern="HistoricalData.csv", full.names=TRUE)
    
      if (length(History == 1)){
      HistoricalData <- read.csv(HistoricalPath, check.names=FALSE)
      HistoricalData$Date <- lubridate::ymd(HistoricalData$Date)
      if (ncol(ShinyData) == ncol(HistoricalData)){
        TheShiniestData <- bind_rows(ShinyData, HistoricalData)
        write.csv(TheShiniestData, HistoricalPath, row.names = FALSE)
      } else {stop("Shiny Historical Data Conflicting Column Numbers")}
      } else {write.csv(ShinyData, HistoricalPath, row.names = FALSE)}
      

    # Regular Order
    TheArchive <- file.path(Folder, "Archive")
    ArchivedDataFile <- list.files(TheArchive, pattern="Archived",
                                   full.names = TRUE)

    if (!length(ArchivedDataFile)==0){

      if(length(ArchivedDataFile)==1){
        ArchivedData <- read.csv(ArchivedDataFile[1], check.names=FALSE)
      } else {message("Two csv files in the folder found!")}

      ArchivedData$DateTime <- lubridate::ymd_hms(ArchivedData$DateTime)
      ArchivedData <- ArchivedData |> mutate(across(starts_with("Flag"), ~ as.logical(.)))

      # Troubleshooting
      if (!ncol(ArchivedData) == ncol(Parsed)){
        Recent <- setdiff(colnames(ArchivedData), colnames(Parsed))
        Previous <- setdiff(colnames(Parsed), colnames(ArchivedData))

        if (length(Previous) == 0){
          UpToHere <- nrow(Parsed)
          WorkAround <- bind_rows(Parsed, ArchivedData)
          WorkAround1 <- WorkAround[1:UpToHere,]
          NewData <- generics::setdiff(WorkAround1, ArchivedData)
          UpdatedData <- rbind(NewData, ArchivedData)
        } else {stop("Mismatched Columns, newer data fewer columns than old data")}
      } else{
        NewData <- generics::setdiff(Parsed, ArchivedData)
        UpdatedData <- rbind(NewData, ArchivedData)
      }

      file.remove(ArchivedDataFile)

    } else {UpdatedData <- Parsed}

    file.remove(DailyQCFiles)

    UpdatedData <- UpdatedData |> arrange(desc(DateTime))

    name <- paste0("ArchivedData", x, ".csv")
    StorageLocation <- file.path(TheArchive, name)
    write.csv(UpdatedData, StorageLocation, row.names=FALSE)
  } else {message("No DailyQCFiles files to update with in ", x)}

}