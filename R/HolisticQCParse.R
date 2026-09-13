#' Dashboard Internal, updates Gain/MFI tracking for Holistic Data CSV
#'
#' @param x The cytometer folder name
#' @param MainFolder The file.path to main folder
#' @param Template Default NULL, a file.path to an openCyto gating template want to apply fist.
#' @param subsets Default NULL, a GatingHierarchy subset to retrieve information from. 
#' @param FuckIt Default FALSE, when the user messed with settings so badly to cause a migraine. 
#' @param sample.name Default '$NAME', used for metadata. 
#'
#' @importFrom flowWorkspace load_cytoset_from_fcs GatingSet keyword
#' @importFrom purrr map
#' @importFrom dplyr bind_rows mutate relocate arrange desc anti_join
#' @importFrom utils read.csv write.csv
#' @importFrom lubridate ymd_hms ymd hms
#' @importFrom openCyto gatingTemplate gt_gating
#' @importFrom data.table fread
#' @importFrom Biobase exprs
#' @importFrom stringr str_c str_sub
#'
#' @return Updated MFI tracking CSV
#' @noRd
HolisticQCParse <- function(x,
                             MainFolder,
                             Template = NULL,
                             subsets = NULL,
                             FuckIt = FALSE,
                             sample.name = "$DATE") {

  Folder <- file.path(MainFolder, x)
  FCS_Files <- list.files(Folder, pattern = "fcs", full.names = TRUE)

  if (!length(FCS_Files) == 0) {

    if (FuckIt == TRUE) {
      Screen <- CytosetScreen(files = FCS_Files)
      MainList <- which.max(sapply(Screen, length))
      Screen <- Screen[MainList][[1]]
      The_CS <- load_cytoset_from_fcs(files = Screen, transformation = FALSE,
                                       truncate_max_range = FALSE)
    } else {
      The_CS <- load_cytoset_from_fcs(files = FCS_Files,
                                       transformation = FALSE,
                                       truncate_max_range = FALSE)
    }

    DateFormat <- keyword(The_CS[[1]])$`$DATE`

    if (sample.name == "$DATE" && DateFormat == "01-Jan-0001") {
      sample.name1 <- "$FIL"
    } else {
      sample.name1 <- sample.name
    }

    if (is.null(Template)) {
      Parsed <- map(.x = The_CS, .f = QC_GainMonitoring,
                    sample.name = sample.name, stats = "median") |>
        bind_rows()
    } else {
      Gating <- fread(Template)
      MyGatingSet <- GatingSet(The_CS)
      MyGatingTemplate <- gatingTemplate(Gating)

      MyGatingSet <- GateCheck(gs = MyGatingSet,
                                gatingtemplate = MyGatingTemplate,
                                subsets = subsets)

      if (is.null(MyGatingSet)) {
        return(MyGatingSet)
      }

      gt_gating(MyGatingTemplate, MyGatingSet)

      Parsed <- map(.x = MyGatingSet, .f = QC_GainMonitoring,
                    subsets = subsets, sample.name = sample.name1,
                    stats = "median") |>
        bind_rows()

      if (DateFormat == "01-Jan-0001") {
        Parsed <- Parsed |>
          mutate(DATE = gsub("DailyQCDataSample_", "", SAMPLE))
        Parsed <- Parsed |> mutate(TIME = sub("^[^_]*_", "", DATE))
        Parsed <- Parsed |> mutate(DATE = sub("_.*", "", DATE))
        Parsed$DATE <- ymd(Parsed$DATE)
        Parsed <- Parsed |>
          mutate(TIME = str_c(str_sub(TIME, 1, 2), ":",
                               str_sub(TIME, 3, 4), ":",
                               str_sub(TIME, 5, 6)),
                 TIME = hms(TIME))
      }

    }

    Parsed <- Parsed |>
      mutate(DateTime = DATE + TIME) |>
      relocate(DateTime, .before = DATE)

    Parsed <- Parsed |> arrange(desc(DateTime))

    ArchiveFolder <- file.path(Folder, "Archive")
    ArchiveCSV <- list.files(ArchiveFolder, pattern = "Holistic",
                              full.names = TRUE)

    if (!length(ArchiveCSV) == 0) {

      if (!length(ArchiveCSV) > 1) {

        ArchiveData <- read.csv(ArchiveCSV, check.names = FALSE)
        ArchiveData$DateTime <- ymd_hms(ArchiveData$DateTime)
        ArchiveData$DATE <- ymd(ArchiveData$DATE)
        ArchiveData$TIME <- hms(ArchiveData$TIME)

        if (!ncol(Parsed) == ncol(ArchiveData)) {
          message("Mismatched Number of Columns")
        }

        NewData <- Parsed |>
          anti_join(ArchiveData, by = c("DATE", "TIME"))

        UpdatedData <- bind_rows(NewData, ArchiveData)

        file.remove(ArchiveCSV)

      } else {
        stop("Two Holistic csv files in the archive folder!")
      }

    } else {
      UpdatedData <- Parsed
    }

    UpdatedData <- UpdatedData |> arrange(desc(DateTime))

    file.remove(FCS_Files)
    name <- paste0("HolisticData", x, ".csv")
    StorageLocation <- file.path(ArchiveFolder, name)
    write.csv(UpdatedData, StorageLocation, row.names = FALSE)

  } else {
    message("No fcs files to update with in ", x)
  }
}