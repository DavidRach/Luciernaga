utils::globalVariables(c("%rCV", ".", ".data", "AdjustedY", "AggregateName",
                         "Area Scaling Factor", "Averaged", "Backups", "Baseline",
                         "Brightness", "CellRatio_Retained", "Cluster", "Clusters",
                          "ColorSelection", "Column", "ConcentrationScientific",
                         "Condition", "Count", "Counts", "Creator",
                         "CurrentConcentration", "DATE", "Date", "DateTime", "Decision",
                         "DeltaGain", "DesiredConcentration", "Detector", "Detector1Raw",
                         "Detectors", "Experiment", "Fluorochrome", "Fluorophore",
                         "Fluors", "Freq", "From", "Gain", "GainFlag", "Gain_Logical",
                          "IncreaseVolumeML", "Instrument", "Laser", "Laser Delay",
                         "Laser Power", "Ligand", "LostRatio", "LowerLimit_MillionCells",
                         "MainDetector", "NewID", "Partition", "Percentiles", "QCStatus",
                         "RCVFlag", "RankImportance", "RankValue", "Ratio", "Sample",
                          "Specimen", "SpinDown", "Splitpoint", "Summed", "TIME",
                         "TheDetector", "TheDetectors", "TheHeight", "TheInstrumentCount",
                         "TheMedian", "TheSampleName", "TheSamples", "TheSums",
                         "Theoretical", "Time", "Timepoint", 'To', "TotalCells",
                         "TotalTubes", "TotalVolume", "Type", "Var1", "Var2",
                         "cell_fill", "cell_text", "cols_label", "cutoff", "dims",
                          "matches", "name", "opt_table_outline", "parent",
                         "px", "rCV", "rCV_Logical", "regular_split", "specimen",
                         "tab_options", "tags", "value", "x", "xlim", "yhat",
                         "ylim", "pData<-", "parameters<-"))

#' Small Internal Function
#' 
#' @param data Something
#' @param x Something
#' @param type Something
#' 
#' @importFrom dplyr select
#' @importFrom tidyselect all_of
#' @importFrom dplyr slice
#' @importFrom dplyr pull
#' 
#' @return Some value
#' 
#' @noRd
CurrentStatus <- function(data, x, type){
  Status <- data |> select(all_of(x)) |> na.omit() |>
    slice(1)|> pull()
  return(Status)
}


#' Small Internal Function
#' 
#' @param x Something
#' 
#' @importFrom dplyr case_when
#' 
#' @return Something
#' 
#' @noRd                       
InstrumentText <- function(x) {
  dplyr::case_when(
    x == "Green" ~ "Pass",
    x == "Yellow" ~ "Caution",
    x == "Orange" ~ "Caution",
    x == "Red" ~ "Fail",
    TRUE ~ NA_character_)
}

#' Small Internal Function
#' 
#' @param x Something
#' 
#' @importFrom dplyr case_when
#' 
#' @return Something
#' 
#' @noRd
InstrumentColor <- function(x) {
  dplyr::case_when(
    x == "Green" ~ "success",
    x == "Yellow" ~ "caution",
    x == "Orange" ~ "warning",
    x == "Red" ~ "danger",
    TRUE ~ NA_character_)
}

#' Dashboard Internal, updates archive ApplicationLog.csv from Setup
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
#'
#' @return Updated tracking data CSV in the Archive Folder
#' @noRd
AppQCParse <- function(MainFolder, x){

  Folder <- file.path(MainFolder, x)
  DailyQCFiles <- list.files(Folder, pattern="Application",
                                full.names = TRUE)

  if (!length(DailyQCFiles)==0){

    if (length(DailyQCFiles)>=1){

      Parsed <- map(.x=DailyQCFiles, .f=Luciernaga:::ApplicationLogParse) %>% bind_rows()
      Parsed <- Parsed[grepl("cmd: SitFlush: Begin", Parsed$Command), ]
      Parsed <- Parsed %>% arrange(desc(DateTime))

    } else {stop("Two csv files in the folder found!")}

    TheArchive <- file.path(Folder, "Archive")
    ArchivedDataFile <- list.files(TheArchive, pattern="Application",
                                   full.names = TRUE)

    if (!length(ArchivedDataFile)==0){

      if(length(ArchivedDataFile)==1){
        ArchivedData <- read.csv(ArchivedDataFile[1], check.names=FALSE)
      } else {message("Two csv files in the folder found!")}

      # Troubleshooting
      if (!ncol(ArchivedData) == ncol(Parsed)){
        stop("Mismatched Number of Columns")
        }

      ArchivedData$DateTime <- lubridate::ymd_hms(ArchivedData$DateTime)
      ArchivedData <- ArchivedData |> arrange(desc(DateTime))
      NewData <- generics::setdiff(Parsed, ArchivedData)
      UpdatedData <- rbind(NewData, ArchivedData)

      file.remove(ArchivedDataFile)

    } else {UpdatedData <- Parsed}

    file.remove(DailyQCFiles)

    UpdatedData <- UpdatedData %>% arrange(desc(DateTime))

    name <- paste0("ApplicationData", x, ".csv")
    StorageLocation <- file.path(TheArchive, name)
    write.csv(UpdatedData, StorageLocation, row.names=FALSE)
  } else {message("No DailyQCFiles files to update with in ", x)}

}

#' Internal for Dashboard, plots User Usage for respective instruments all time
#' 
#' @param data The Data derrived from AppParse filtered for SitFlushes
#' @param TheInstrument Instrument designation in the Instrument column, example "5L"
#' 
#' @importFrom dplyr filter
#' @importFrom lubridate wday
#' @importFrom lubridate hour
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom dplyr arrange
#' @importFrom dplyr desc
#' @importFrom dplyr pull
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 geom_col
#' @importFrom ggplot2 facet_grid
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 scale_x_continuous
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_text
#' 
#' @return A ggplot2 object
#' 
#' @noRd
UsagePlot <- function(data, TheInstrument){
  data <- data |> filter(Instrument %in% TheInstrument)
  
  data$DayOfWeek <- wday(data$DateTime, label = TRUE, abbr = TRUE)
  data$DayOfWeek <- factor(data$DayOfWeek,
   levels = c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun"))

  data$Hour <- hour(data$DateTime)

  dataByHourDay <- data |>
  group_by(DayOfWeek, Hour) |>
  summarise(count = n(), .groups = "drop")
  
  MaxCount <- dataByHourDay |> arrange(desc(count)) |>
    slice(1) |> pull(count)
  MaxCount <- MaxCount*1.05
  MaxCount <- round(MaxCount, 0)

  TheTitle <- paste0("Aurora ", TheInstrument, " All Time Usage")

  plot <- ggplot(dataByHourDay, aes(x = Hour, y = count)) +
    geom_col(fill = "black") + 
    labs(title = TheTitle,
        x = "Hour of Day",
        y = "Sit Flushes") +
    facet_grid(DayOfWeek ~ ., scales = "free_y") +
    theme_bw() +
    scale_x_continuous(breaks = 0:23) + lims(y=c(NA, MaxCount)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text.y = element_text(angle = 0)) 

  return(plot)
}

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
      HistoricalData <- read.csv(HistoricalPath, check.names=FALSE)
      HistoricalData$Date <- lubridate::ymd(HistoricalData$Date)
      if (ncol(ShinyData) == ncol(HistoricalData)){
        TheShiniestData <- bind_rows(ShinyData, HistoricalData)
        write.csv(TheShiniestData, HistoricalPath, row.names = FALSE)
      } else {stop("Shiny Historical Data Conflicting Column Numbers")}

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

#' Dashboard Internal, updates instrument tracking CSV
#'
#' @param MainFolder The file.path to the Main Folder
#' @param x The Cytometer Folder Name
#' @param Maintainer Logical override for when number of columns don't match
#'
#' @importFrom dplyr mutate
#' @importFrom dplyr across
#' @importFrom tidyselect starts_with
#' @importFrom utils read.csv
#' @importFrom lubridate ymd_hms
#' @importFrom generics setdiff
#' @importFrom utils write.csv
#'
#' @return Updated tracking data CSV in the Archive Folder
#' @noRd
LevyJenningsParse <- function(MainFolder, x, Maintainer=FALSE){

  Folder <- file.path(MainFolder, x)
  LJTrackingFiles <- list.files(Folder, pattern="LJ",
                                full.names = TRUE)

  if (!length(LJTrackingFiles)==0){

    if (length(LJTrackingFiles)==1){

      Parsed <- QC_FilePrep(x=LJTrackingFiles, TrackChange=FALSE)
      Parsed <- Parsed %>% mutate(across(starts_with("Flag"), ~ as.logical(.)))

    } else {stop("Two csv files in the folder found!")}

    TheArchive <- file.path(Folder, "Archive")
    ArchivedDataFile <- list.files(TheArchive, pattern="Archived",
                                   full.names = TRUE)

    if (!length(ArchivedDataFile)==0){

      if(length(ArchivedDataFile)==1){
        ArchivedData <- read.csv(ArchivedDataFile, check.names=FALSE)
      } else {message("Two csv files in the folder found!")}

      # Troubleshooting
      if (!ncol(ArchivedData) == ncol(Parsed)){

        if (Maintainer==TRUE){
          TheseColumns <- setdiff(colnames(Parsed), colnames(ArchivedData))

          for (col in TheseColumns) {
            ArchivedData[[col]] <- NA
          }

          if (!ncol(ArchivedData) == ncol(Parsed)){stop("Still no rescue")}

        } else {
          stop("The number of columns for the new data don't match
       that of the archived data. Please make sure to
       export the Levy-Jennings trackings with all available
       parameters")
        }
      }

      ArchivedData$DateTime <- ymd_hms(ArchivedData$DateTime)
      ArchivedData <- ArchivedData %>% mutate(across(starts_with("Flag"), ~ as.logical(.)))
      NewData <- generics::setdiff(Parsed, ArchivedData)
      UpdatedData <- rbind(NewData, ArchivedData)

      file.remove(ArchivedDataFile)

    } else {UpdatedData <- Parsed}

    file.remove(LJTrackingFiles)

    name <- paste0("ArchivedData", x, ".csv")
    StorageLocation <- file.path(TheArchive, name)
    write.csv(UpdatedData, StorageLocation, row.names=FALSE)
  } else {message("No LevyJennings files to update with in ", x)}

}

#' Dashboard Internal, updates MFI tracking CSV
#'
#' @param x The cytometer folder name
#' @param MainFolder The file.path to main folder
#'
#' @importFrom flowWorkspace load_cytoset_from_fcs
#' @importFrom purrr map
#' @importFrom dplyr bind_rows
#' @importFrom dplyr mutate
#' @importFrom dplyr relocate
#' @importFrom dplyr arrange
#' @importFrom dplyr desc
#' @importFrom utils read.csv
#' @importFrom lubridate ymd_hms
#' @importFrom lubridate ymd
#' @importFrom lubridate hms
#' @importFrom dplyr anti_join
#' @importFrom utils write.csv
#'
#' @return Updated MFI tracking CSV
#' @noRd
QCBeadParse <- function(x, MainFolder){
  Folder <- file.path(MainFolder, x)
  FCS_Files <- list.files(Folder, pattern="fcs", full.names=TRUE)

  if(!length(FCS_Files) == 0){

    QCBeads <- FCS_Files[grep("Before|After", FCS_Files)]

    if (length(QCBeads) == 0){
      message("No Before After detected in names")
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

    BeforeAfter <- map(.x=BeforeAfter_CS, .f=QC_GainMonitoring,
                       sample.name = "TUBENAME", stats="median") %>% bind_rows()

    BeforeAfter <- BeforeAfter %>% mutate(DateTime = DATE+TIME) %>%
      relocate(DateTime, .before=DATE)

    BeforeAfter <- BeforeAfter %>% arrange(desc(DateTime))

    ArchiveFolder <- file.path(Folder, "Archive")
    ArchiveCSV <- list.files(ArchiveFolder, pattern="Bead", full.names=TRUE)

    if (!length(ArchiveCSV) == 0){

    if (!length(ArchiveCSV) > 1){

      ArchiveData <- read.csv(ArchiveCSV, check.names=FALSE)
      ArchiveData$DateTime <- lubridate::ymd_hms(ArchiveData$DateTime)
      ArchiveData$DATE <- lubridate::ymd(ArchiveData$DATE)
      ArchiveData$TIME <- lubridate::hms(ArchiveData$TIME)

      if (!ncol(BeforeAfter) == ncol(ArchiveData)){
        stop("Mismatched Number of Columns")
      }

      NewData <- BeforeAfter %>%
        anti_join(ArchiveData, by = c("DATE", "TIME"))

      UpdatedData <- rbind(NewData, ArchiveData)

      file.remove(ArchiveCSV)

    } else {stop("Two BeadData csv files in the archive folder!")}

    } else {UpdatedData <- BeforeAfter}

    UpdatedData <- UpdatedData %>% arrange(desc(DateTime))

    file.remove(FCS_Files)
    name <- paste0("BeadData", x, ".csv")
    StorageLocation <- file.path(ArchiveFolder, name)
    write.csv(UpdatedData, StorageLocation, row.names=FALSE)

    } else {message("No fcs files to update with in ", x)}
}

#' Dashboard Internal, updates Gain/MFI tracking for Holistic Data CSV
#'
#' @param x The cytometer folder name
#' @param MainFolder The file.path to main folder
#' @param Template Default NULL, a file.path to an openCyto gating template want to apply fist.
#' @param subsets Default NULL, a GatingHierarchy subset to retrieve information from. 
#' @param FuckIt Default FALSE, when the user messed with settings so badly to cause a migraine. 
#'
#' @importFrom flowWorkspace load_cytoset_from_fcs
#' @importFrom purrr map
#' @importFrom dplyr bind_rows
#' @importFrom dplyr mutate
#' @importFrom dplyr relocate
#' @importFrom dplyr arrange
#' @importFrom dplyr desc
#' @importFrom utils read.csv
#' @importFrom lubridate ymd_hms
#' @importFrom lubridate ymd
#' @importFrom lubridate hms
#' @importFrom dplyr anti_join
#' @importFrom utils write.csv
#' @importFrom flowWorkspace GatingSet
#' @importFrom openCyto gatingTemplate
#' @importFrom openCyto gt_gating
#' @importFrom data.table fread
#' @importFrom Biobase exprs
#'
#' @return Updated MFI tracking CSV
#' @noRd
HolisticQCParse <- function(x, MainFolder, Template=NULL, subsets=NULL, FuckIt=FALSE){
  Folder <- file.path(MainFolder, x)
  FCS_Files <- list.files(Folder, pattern="fcs", full.names=TRUE)

  if(!length(FCS_Files) == 0){

    if (FuckIt == TRUE){
    Screen <- CytosetScreen(files = FCS_Files)
    MainList <- which.max(sapply(Screen, length))
    Screen <- Screen[MainList][[1]]
    The_CS <- load_cytoset_from_fcs(files = Screen,
       transformation = FALSE, truncate_max_range = FALSE)
    } else {
      The_CS <- load_cytoset_from_fcs(
        files = FCS_Files, transformation = FALSE,
        truncate_max_range = FALSE)
    }

    if (is.null(Template)){
    Parsed <- map(.x=The_CS, .f=QC_GainMonitoring,
                       sample.name = "$DATE", stats="median") |> bind_rows()
    } else {
      Gating <- data.table::fread(Template)
      MyGatingSet <- GatingSet(The_CS)
      MyGatingTemplate <- gatingTemplate(Gating)

      MyGatingSet <- GateCheck(gs=MyGatingSet, gatingtemplate = MyGatingTemplate)
      
      if (is.null(MyGatingSet)){return(MyGatingSet)}

      gt_gating(MyGatingTemplate, MyGatingSet)

    Parsed <- map(.x=MyGatingSet, .f=QC_GainMonitoring, subsets=subsets,
        sample.name = "$DATE", stats="median") |> bind_rows()
      
    }

    Parsed <- Parsed |> mutate(DateTime = DATE+TIME) |>
      relocate(DateTime, .before=DATE)

    Parsed <- Parsed |> arrange(desc(DateTime))

    ArchiveFolder <- file.path(Folder, "Archive")
    ArchiveCSV <- list.files(ArchiveFolder, pattern="Holistic", full.names=TRUE)

    if (!length(ArchiveCSV) == 0){

    if (!length(ArchiveCSV) > 1){

      ArchiveData <- read.csv(ArchiveCSV, check.names=FALSE)
      ArchiveData$DateTime <- lubridate::ymd_hms(ArchiveData$DateTime)
      ArchiveData$DATE <- lubridate::ymd(ArchiveData$DATE)
      ArchiveData$TIME <- lubridate::hms(ArchiveData$TIME)

      if (!ncol(Parsed) == ncol(ArchiveData)){
        message("Mismatched Number of Columns")
      }

      NewData <- Parsed %>%
        anti_join(ArchiveData, by = c("DATE", "TIME"))

      UpdatedData <- bind_rows(NewData, ArchiveData)

      file.remove(ArchiveCSV)

    } else {stop("Two Holistic csv files in the archive folder!")}

    } else {UpdatedData <- Parsed}

    UpdatedData <- UpdatedData %>% arrange(desc(DateTime))

    file.remove(FCS_Files)
    name <- paste0("HolisticData", x, ".csv")
    StorageLocation <- file.path(ArchiveFolder, name)
    write.csv(UpdatedData, StorageLocation, row.names=FALSE)

  } else {message("No fcs files to update with in ", x)}
}

#' Similar to CytosetScreen, checks for mismatching cytoframes that throw inconvenient errors
#' 
#' @param gs A gating set object
#' @param gatingtemplate The gating template object
#' 
#' @return Purified Gating Set or a NULL Value
#' 
#' @noRd
GateCheck <- function(gs, gatingtemplate){
  
  Nodes <- gatingtemplate@nodes
  FinalNode <- Nodes[length(Nodes)]
  These <- Nodes[-1]
  CheckThis <- paste(These, collapse="|")
  ThisFluor <- gatingtemplate@edgeData@data[[CheckThis]]$gtMethod@dims

  TheIndex <- any(colnames(gs) %in% ThisFluor)
  Present <- gs[TheIndex]

  if (length(Present) == 0){gs <- NULL
  } else {gs <- Present}
  return(gs)
}

#' Dashboard Internal, updates MFI tracking CSV
#'
#' @param x The cytometer folder name
#' @param MainFolder The file.path to main folder
#'
#' @importFrom flowWorkspace load_cytoset_from_fcs
#' @importFrom purrr map
#' @importFrom dplyr bind_rows
#' @importFrom dplyr mutate
#' @importFrom dplyr relocate
#' @importFrom dplyr arrange
#' @importFrom dplyr desc
#' @importFrom utils read.csv
#' @importFrom lubridate ymd_hms
#' @importFrom lubridate ymd
#' @importFrom lubridate hms
#' @importFrom dplyr anti_join
#' @importFrom utils write.csv
#'
#' @return Updated MFI tracking CSV
#' @noRd
QCBeadParse <- function(x, MainFolder){
  Folder <- file.path(MainFolder, x)
  FCS_Files <- list.files(Folder, pattern="fcs", full.names=TRUE)

  if(!length(FCS_Files) == 0){

    QCBeads <- FCS_Files[grep("Before|After", FCS_Files)]
    BeforeAfter_CS <- load_cytoset_from_fcs(files=QCBeads,
                                            transformation=FALSE, truncate_max_range = FALSE)

    BeforeAfter <- map(.x=BeforeAfter_CS, .f=QC_GainMonitoring,
                       sample.name = "TUBENAME", stats="median") %>% bind_rows()

    BeforeAfter <- BeforeAfter %>% mutate(DateTime = DATE+TIME) %>%
      relocate(DateTime, .before=DATE)

    BeforeAfter <- BeforeAfter %>% arrange(desc(DateTime))

    ArchiveFolder <- file.path(Folder, "Archive")
    ArchiveCSV <- list.files(ArchiveFolder, pattern="Bead", full.names=TRUE)

    if (!length(ArchiveCSV) == 0){

    if (!length(ArchiveCSV) > 1){

      ArchiveData <- read.csv(ArchiveCSV, check.names=FALSE)
      ArchiveData$DateTime <- lubridate::ymd_hms(ArchiveData$DateTime)
      ArchiveData$DATE <- lubridate::ymd(ArchiveData$DATE)
      ArchiveData$TIME <- lubridate::hms(ArchiveData$TIME)

      if (!ncol(BeforeAfter) == ncol(ArchiveData)){
        stop("Mismatched Number of Columns")
      }

      NewData <- BeforeAfter %>%
        anti_join(ArchiveData, by = c("DATE", "TIME"))

      UpdatedData <- rbind(NewData, ArchiveData)

      file.remove(ArchiveCSV)

    } else {stop("Two BeadData csv files in the archive folder!")}

    } else {UpdatedData <- BeforeAfter}

    UpdatedData <- UpdatedData %>% arrange(desc(DateTime))

    file.remove(FCS_Files)
    name <- paste0("BeadData", x, ".csv")
    StorageLocation <- file.path(ArchiveFolder, name)
    write.csv(UpdatedData, StorageLocation, row.names=FALSE)

    } else {message("No fcs files to update with in ", x)}
}

#' Dashboard Internal, loads updated data
#'
#' @param x The Cytometer Folder Name
#' @param MainFolder The file.path to the main folder
#' @param type Whether to return "MFI" or "Gain" plots
#'
#' @importFrom utils read.csv
#' @importFrom lubridate ymd_hms
#' @importFrom lubridate ymd
#' @importFrom lubridate hms
#' @importFrom stringr str_detect
#' @importFrom lubridate mdy_hm
#' @importFrom dplyr arrange
#' @importFrom dplyr desc
#'
#' @return Updated Tracking Data CSV for specified type
#' @noRd
CurrentData <- function(x, MainFolder, type){

  ArchiveLocation <- file.path(MainFolder, x, "Archive")

  if (type == "MFI"){
    BeadData <- list.files(ArchiveLocation, pattern="Bead",
                           full.names=TRUE)
    Data <- read.csv(BeadData, check.names=FALSE)
    Data$DateTime <- lubridate::ymd_hms(Data$DateTime)
    Data$DATE <- lubridate::ymd(Data$DATE)
    Data$TIME <- lubridate::hms(Data$TIME)
  }

  if (type == "Gain"){
    ArchiveData <- list.files(ArchiveLocation, pattern="Archived",
                              full.names=TRUE)
    Data <- read.csv(ArchiveData, check.names=FALSE)
    #lubridate::ymd_hms(Data$DateTime)

    if (any(str_detect(Data$DateTime, ":.*:"))){
      Data$DateTime <- lubridate::ymd_hms(Data$DateTime)
    } else {Data$DateTime <- lubridate::mdy_hm(Data$DateTime)}
  }

  if (type == "Both"){
    BothData <- list.files(ArchiveLocation, pattern="Holistic",
                           full.names=TRUE)
    Data <- read.csv(BothData, check.names=FALSE)
    Data$DateTime <- lubridate::ymd_hms(Data$DateTime)
    Data$DATE <- lubridate::ymd(Data$DATE)
    Data$TIME <- lubridate::hms(Data$TIME)
  }

  Data <- Data %>% arrange(desc(DateTime))
  return(Data)
}

#' Dashboard Internal, rearranges vector by color
#'
#' @param colors A vector of cytometer parameters to be rearranged
#'
#' @return Reordered vector according to light wavelength
#' @noRd
ColorPriority <- function(colors){
  Ordered <- colors[order(grepl("^Ultra", colors) * -1,
                          grepl("^UV", colors) * -1,
                          grepl("^Violet", colors) * -1,
                          grepl("^Blue", colors) * -1,
                          grepl("^Yellow", colors) * -1,
                          grepl("^yellow", colors) * -1,
                          grepl("^Red", colors) * -1)]
  return(Ordered)
}

#' Dashboard Internal, rearranges vector by scatter parameter
#'
#' @param colors A vector of scatter parameters to be rearranged
#'
#' @return Rearranged vector of scatter parameters
#' @noRd
ScalePriority <- function(colors){
  Ordered <- colors[order(grepl("^FSC", colors) * -1,
                          grepl("^SSC", colors) * -1,
                          grepl("-A", colors) * -1,
                          grepl("-H", colors) * -1,
                          grepl("-W", colors) * -1)]
  return(Ordered)
}

#' Dashboard Internal, processes to did parameter pass in past week
#'
#' @param x The data.frame output from LevyJennings or QCBeads Parse
#'
#' @importFrom dplyr mutate
#' @importFrom dplyr relocate
#' @importFrom dplyr pull
#' @importFrom purrr map
#' @importFrom dplyr bind_rows
#'
#' @return Data frame of passing status for respective parameters
#' @noRd
ShinyQCSummary <- function(x, Instrument){

  Intermediate <- x %>% mutate(Date=as.Date(DateTime)) %>%
    relocate(Date, .before=1)

  Dates <- Intermediate %>% pull(Date) %>% unique()

  Data <- map(.x=Dates, .f=ShinyQCSummaryParser, Intermediate=Intermediate) %>%
    bind_rows()

  Data <- Data %>% mutate(Instrument=Instrument) %>%
    relocate(Instrument, .after="Date")

  return(Data)
}

#' Dashboard Internal, processes to did parameter pass in past week
#'
#' @param x The iterated date
#' @param Intermediate The original instrument data
#'
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom lubridate weeks
#' @importFrom tidyselect starts_with
#' @importFrom tidyselect contains
#' @importFrom tidyselect all_of
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr left_join
#' @importFrom dplyr pull
#' @importFrom purrr map
#' @importFrom dplyr bind_rows
#'
#' @return Data frame of passing status for respective parameters
#' @noRd
ShinyQCSummaryParser <- function(x, Intermediate){
  Date <- x

  x <- Intermediate %>% filter(Date %in% x)

  if (nrow(x) > 1){
    #Data <- x %>% filter(DateTime > WindowOfInterest)
    Data <- x %>% slice(1)
  } else {Data <- x}

  Flags <- Data %>% select(starts_with("Flag"))
  colnames(Flags) <- gsub("Flag-", "", colnames(Flags))
  Gains <- Flags %>% select(contains("Gain"))
  TheGains <- colnames(Gains)
  rCV <- Flags %>% select(contains("rCV"))
  TherCV <- colnames(rCV)

  TheGainData <- Data %>% select(all_of(c("DateTime", TheGains)))
  colnames(Gains) <- gsub("-Gain", "", fixed=TRUE, colnames(Gains))
  colnames(TheGainData) <- gsub("-Gain", "", fixed=TRUE, colnames(TheGainData))

  TheGainData <- TheGainData %>%
    pivot_longer(!DateTime, names_to = "Detector", values_to = "Gain")

  Gains <- Gains %>% mutate(DateTime=Data$DateTime) %>% relocate(DateTime, .before=1)

  Gains <- Gains %>%
    pivot_longer(!DateTime, names_to = "Detector", values_to = "Gain_Logical")

  TherCVData <- Data %>% select(all_of(c("DateTime", TherCV)))
  colnames(rCV) <- gsub("-% rCV", "", fixed=TRUE, colnames(rCV))
  colnames(TherCVData) <- gsub("-% rCV", "", fixed=TRUE, colnames(TherCVData))

  TherCVData <- TherCVData %>% pivot_longer(!DateTime, names_to = "Detector", values_to = "rCV")

  rCV <- rCV %>% mutate(DateTime=Data$DateTime) %>% relocate(DateTime, .before=1)

  rCV <- rCV %>% pivot_longer(!DateTime, names_to = "Detector", values_to = "rCV_Logical")

  Tidy <- TheGainData %>%
    left_join(Gains, by = c("Detector", "DateTime")) %>%
    left_join(TherCVData, by = c("Detector", "DateTime")) %>%
    left_join(rCV, by = c("Detector", "DateTime"))

  TheDetectors <- Tidy %>% pull(Detector) %>% unique()

  Summary <- map(.x=TheDetectors, .f=QCSummaryCheck, data=Tidy) %>% bind_rows()

  Summary <- Summary %>% mutate(Date=Date) %>% relocate(Date, .before=1)
  return(Summary)
}


#' Dashboard Internal, processes to did parameter pass in past week
#'
#' @param x The data.frame output from LevyJennings or QCBeads Parse
#'
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom lubridate weeks
#' @importFrom tidyselect starts_with
#' @importFrom tidyselect contains
#' @importFrom tidyselect all_of
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr left_join
#' @importFrom dplyr pull
#' @importFrom purrr map
#' @importFrom dplyr bind_rows
#'
#' @return Data frame of passing status for respective parameters
#' @noRd
VisualQCSummary <- function(x){

  WindowOfInterest <- Sys.time() - weeks(1)

  if (nrow(x) > 1){
  Data <- x |> filter(DateTime > WindowOfInterest)

  if(nrow(Data) == 0){Data <- x |> slice(1)}

  } else {Data <- x}

  Flags <- Data |> select(starts_with("Flag"))
  colnames(Flags) <- gsub("Flag-", "", colnames(Flags))
  Gains <- Flags |> select(contains("Gain"))
  TheGains <- colnames(Gains)
  rCV <- Flags |> select(contains("rCV"))
  TherCV <- colnames(rCV)

  TheGainData <- Data |> select(all_of(c("DateTime", TheGains)))
  colnames(Gains) <- gsub("-Gain", "", fixed=TRUE, colnames(Gains))
  colnames(TheGainData) <- gsub("-Gain", "", fixed=TRUE, colnames(TheGainData))

  TheGainData <- TheGainData |>
    pivot_longer(!DateTime, names_to = "Detector", values_to = "Gain")

  Gains <- Gains |> mutate(DateTime=Data$DateTime) |> relocate(DateTime, .before=1)

  Gains <- Gains |>
    pivot_longer(!DateTime, names_to = "Detector", values_to = "Gain_Logical")

  TherCVData <- Data |> select(all_of(c("DateTime", TherCV)))
  colnames(rCV) <- gsub("-% rCV", "", fixed=TRUE, colnames(rCV))
  colnames(TherCVData) <- gsub("-% rCV", "", fixed=TRUE, colnames(TherCVData))

  TherCVData <- TherCVData |> pivot_longer(!DateTime, names_to = "Detector", values_to = "rCV")

  rCV <- rCV |> mutate(DateTime=Data$DateTime) |> relocate(DateTime, .before=1)

  rCV <- rCV |> pivot_longer(!DateTime, names_to = "Detector", values_to = "rCV_Logical")

  Tidy <- TheGainData |>
    left_join(Gains, by = c("Detector", "DateTime")) |>
    left_join(TherCVData, by = c("Detector", "DateTime")) |>
    left_join(rCV, by = c("Detector", "DateTime"))

  TheDetectors <- Tidy |> pull(Detector) |> unique()

  Summary <- map(.x=TheDetectors, .f=QCSummaryCheck, data=Tidy) |> bind_rows()
  return(Summary)
}

#' Dashboard Internal, summarizes passing status
#'
#' @param x The iterated parameter
#' @param data The data being compared against
#'
#' @importFrom dplyr filter
#' @importFrom dplyr slice
#' @importFrom dplyr pull
#'
#' @return The color-coded summary
#' @noRd
QCSummaryCheck <- function(x, data){
  Subset <- data |> dplyr::filter(Detector %in% x)

  GainValue <- Subset |> slice(1) |> pull(Gain)

  if (any(Subset$Gain_Logical == TRUE)){
    Followup <- Subset |> slice(1) |> pull(Gain_Logical)
    if (Followup == TRUE){GainStatus <- "Red"
    } else {GainStatus <- "Yellow"}
  } else {GainStatus <- "Green"}

  rCVValue <- Subset |> slice(1) |> pull(rCV) |> round(2)

  if (any(Subset$rCV_Logical == TRUE)){
    Followup <- Subset |> slice(1) |> pull(rCV_Logical)
    if (Followup == TRUE){rCVStatus <- "Red"
    } else {rCVStatus <- "Yellow"}
  } else {rCVStatus <- "Green"}

  Summary <- data.frame(Detector=x, GainValue=GainValue, Gain=GainStatus,
    rCVValue = rCVValue, rCV=rCVStatus)
  return(Summary)
}

#' Dashboard Internal, holistic passing summary
#'
#' @param x The Instrument Name
#' @param y The Instrument Data
#'
#' @return Global Passing Status
#' @noRd
ColorCodeStatus <- function(x, y){

  data <- y

  CytekMandate <- c("FSC", "SSC", "SSC-B")
  CytekData <- data %>%
    filter(Detector %in% CytekMandate)
  RCVdata <- data %>%
    filter(str_detect(Detector, "3")) %>%
    filter(!str_detect(Detector, "1"))
  RCVdata <- rbind(CytekData, RCVdata)

  if(nrow(data)== 0){ColorCode <- "Unknown"}

  if (nrow(data) > 0){

  if (any(data$Gain == "Red") || any(RCVdata$rCV == "Red")) {
    ColorCode <- "Red" # Overall QC Fail
  } else if (any(data$rCV == "Red")) {
    ColorCode <- "Orange" # Non-primary RCV Fail
  } else if (any(data$Gain == "Yellow") || any(data$rCV == "Yellow")) {
    ColorCode <- "Yellow"
  } else {ColorCode <- "Green"}

  }

  QCResults <- data.frame(Instrument = x, QCStatus=ColorCode)
  return(QCResults)
}

#' Dashboard Internal, returns designated hex color.
#'
#' @param x The instrument designation
#' @param data The QC status data.frame
#'
#' @importFrom dplyr filter
#' @importFrom dplyr pull
#'
#' @return A hex code to fill with
#' @noRd
ColorCode <- function(x, data){

  Color <- data %>% dplyr::filter(Instrument %in% x) %>%
    pull(QCStatus)

  Hex <- "#FFFFFF"

  if (Color == "Red"){Hex <- "#C80815"}
  if (Color == "Orange"){Hex <- "#FF6E00"}
  if (Color == "Yellow"){Hex <- "#BA8E23"}
  if (Color == "Green"){Hex <- "#0B6623"}
  return(Hex)
}

#' Dashboard Internal, returns gt table
#'
#' @param data The QC color status returns
#'
#' @importFrom gt gt
#' @importFrom gt data_color
#' @importFrom scales col_factor
#' @importFrom gt sub_values
#' @importFrom gt opt_table_font
#' @importFrom gt cols_align
#' @importFrom gt tab_spanner
#' @importFrom gt cols_label
#' @noRd
SmallTable <- function(data){

  table <- data %>%
    gt() %>%
    data_color(
      columns = c(Gain, rCV),
      fn=function(x) {
        dplyr::case_when(
          x == "Green" ~ "#0B6623",
          x == "Orange" ~ "#FF6E00",
          x == "Yellow" ~ "#BA8E23",
          x == "Red" ~ "#C80815",
          TRUE ~ NA_character_
        )
      }
  )

    Substituted <- table  |>
      sub_values(values= c("Green"), replacement = "Pass") |>
      sub_values(values= c("Orange"), replacement = "Warning") |>
      sub_values(values= c("Yellow"), replacement = "Caution") |>
      sub_values(values= c("Red"), replacement = "Fail")

    Bolded <- Substituted |>
      opt_table_font(font = "Montserrat") |>
      cols_align(align = "center")

    Final <- Bolded |> tab_spanner(
      label = "Gain ",
      columns = c(GainValue, Gain)
    ) |> tab_spanner(
      label = "%RCV ",
      columns = c(rCVValue, rCV)
    ) |> cols_label(
      GainValue = "Value",
      Gain = "Status"
    ) |> cols_label(
      rCVValue = "Value",
      rCV = "Status"
    )

  return(Final)
}

#' Dashboard Internal, summarizes 3 months QC data for all instruments
#'
#' @param x A vector of instrument names
#' @param y A list of LevyJenningParse updated data objects
#' @param timewindow The number  desired months
#' 
#' @importFrom purrr map2
#' @importFrom dplyr bind_rows
#' @importFrom dplyr mutate
#' @importFrom dplyr filter
#' @importFrom dplyr pull
#' @importFrom dplyr group_by
#' @importFrom dplyr slice
#' @importFrom dplyr ungroup
#' @importFrom tidyr pivot_wider
#'
#' @return Data ready for gt coloring
#' @noRd
QCHistory <- function(x, y, timewindow=24){
  InstrumentLength <- length(y)

  if (InstrumentLength > 1){TheInstrumentLength <- 2
  } else {TheInstrumentLength <- 1}

  TheDataset <- map2(.x=x, .f=AcrossTime, .y=y) |> bind_rows()

  TheDates <- TheDataset |> group_by(DateTime) |>
    mutate(TheInstrumentCount = n()) |>
    filter(TheInstrumentCount >= TheInstrumentLength) |> pull(DateTime)

  Assembled <- TheDataset |> filter(DateTime %in% TheDates)

  Assembled <- Assembled |> group_by(DateTime, Instrument) |> slice(1) |> ungroup()

  Figure <- Assembled |> group_by(Instrument) |>
    pivot_wider(names_from = DateTime, values_from = QCStatus)

  return(Figure)
}

#' Dashboard Internal, summarizes x months QC data for all instruments
#'
#' @param x A vector of instrument names
#' @param y A list of LevyJenningParse updated data objects
#' @param timewindow The number  desired months
#' 
#' @importFrom purrr map
#' @importFrom dplyr bind_rows
#' @importFrom dplyr group_by
#' @importFrom dplyr mutate
#' @importFrom dplyr filter
#' @importFrom dplyr pull
#' @importFrom dplyr slice
#' @importFrom dplyr ungroup
#' @importFrom tidyr pivot_wider
#'
#' @return Data ready for gt coloring
#' @noRd
QCHistoryArchive <- function(x, historydata, timewindow=24){
  InstrumentLength <- length(x)

  if (InstrumentLength > 1){TheInstrumentLength <- 2
  } else {TheInstrumentLength <- 2}

  TheDataset <- map(.x=x, data=historydata,
     .f=InternalColorCodeStatus) |> bind_rows()

  TheDates <- TheDataset |> group_by(DateTime) |>
    mutate(TheInstrumentCount = n()) |>
    filter(TheInstrumentCount >= TheInstrumentLength) |> pull(DateTime)

  Assembled <- TheDataset |> filter(DateTime %in% TheDates)

  Assembled <- Assembled |> group_by(DateTime, Instrument) |> slice(1) |> ungroup()

  Figure <- Assembled |> group_by(Instrument) |>
    pivot_wider(names_from = DateTime, values_from = QCStatus)

  return(Figure)
}

#' Internal for QCHistoryArchive
#' 
#' @param x Something
#' @param TheInstrument Something
#' @param TheSubset Something
#' 
#' @importFrom dplyr filter
#' @importFrom dplyr pull
#' @importFrom purrr map
#' @importFrom dplyr bind_rows
#' 
#' @return A value of some form
#' 
#' @noRd
InternalColorCodeStatus <- function(x, data){
  TheInstrument <- x
  TheSubset <- data |> filter(Instrument %in% TheInstrument)
  TheDates <- data |> pull(Date) |> unique()

  TheInstrumentHistory <- map(.x=TheDates, .f=InternalColorDateFilter,
     TheInstrument=TheInstrument, TheSubset=TheSubset) |> bind_rows()
  return(TheInstrumentHistory)
}

#' Internal for QCHistoryArchive
#' 
#' @param x Something
#' @param TheInstrument Something
#' @param TheSubset Something
#' 
#' @importFrom dplyr filter
#' 
#' @return A value of some form
#' 
#' @noRd
InternalColorDateFilter <- function(x, TheInstrument, TheSubset){
  DateTime <- as.Date(x)
  TheDateSummary <- TheSubset |> filter(Date %in% DateTime)
  InstrumentStatus <- ColorCodeStatus(x=TheInstrument, y=TheDateSummary)
  Snapshot <- cbind(DateTime, InstrumentStatus)
  Snapshot$DateTime <- as.Date(Snapshot$DateTime)
  return(Snapshot)
}

#' Dashboard Internal, wrapper for individual instrument history
#'
#' @param x  The Instrument name
#' @param y The Instrument data
#' @param timewindow The number  desired months
#'
#' @importFrom dplyr filter
#' @importFrom dplyr pull
#' @importFrom purrr map
#' @importFrom dplyr bind_rows
#'
#' @return Individual instrument QC history summary
#' @noRd
AcrossTime <- function(x, y, timewindow){
  WindowOfInterest <- Sys.time() - months(timewindow)
  data <- y
  data <- data |> filter(DateTime >= WindowOfInterest)
  TheDates <- data |> pull(DateTime) |> unique()

  Instrument <- x

  # x <- TheDates[1]

  InstrumentHistory <- map(.x=TheDates, data=data, .f=DateMapper, Instrument=Instrument) |> bind_rows()
  return(InstrumentHistory)
}

#' Dashboard Mapper, individual date history retrieval
#'
#' @param x The individual date
#' @param data The instrument data
#' @param Instrument The name of the instrument
#'
#' @importFrom dplyr filter
#'
#' @return Individual date QC summary
#' @noRd
DateMapper <- function(x, data, Instrument){
  DateTime <- x
  TheDateData <- data |> dplyr::filter(DateTime %in% x)
  TheDateSummary <- VisualQCSummary(x=TheDateData)
  InstrumentStatus <- ColorCodeStatus(x=Instrument, y=TheDateSummary)
  Snapshot <- cbind(DateTime, InstrumentStatus)
  Snapshot$DateTime <- as.Date(Snapshot$DateTime)
  return(Snapshot)
}

#' Dashboard Internal, fills color code for global view
#'
#' @param data Assembled data for cytometer pass fails past three months
#'
#' @importFrom gt gt
#' @importFrom gt data_color
#' @importFrom tidyselect everything
#' @importFrom scales col_factor
#' @importFrom gt sub_values
#' @importFrom gt opt_table_font
#' @importFrom gt cols_align
#' @noRd
SmallTableGlobal <- function(data){
  table <- data %>%
    gt() %>%
    data_color(
      columns = c(everything()),
      fn=function(x) {
        x <- as.character(x)

        dplyr::case_when(
          x == "Green" ~ "#0B6623",
          x == "Orange" ~ "#BA8E23",
          x == "Yellow" ~ "#BA8E23",
          x == "Red" ~ "#C80815",
          is.na(x) ~ "#ECECEC",
          TRUE ~ "#FFFFFF"
        )
      }
    )

  Substituted <- table  |>
    sub_values(values= c("Green"), replacement = "Pass") |>
    sub_values(values= c("Orange"), replacement = "Caution") |>
    sub_values(values= c("Yellow"), replacement = "Caution") |>
    sub_values(values= c("Red"), replacement = "Fail")

  Bolded <- Substituted |>
    opt_table_font(font = "Montserrat") |>
    cols_align(align = "center")

  Final <- Bolded

  return(Final)
}
