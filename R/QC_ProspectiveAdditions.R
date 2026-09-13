#' From existing panel, figures out open detectors, and returns potential
#'  fluorophores that might fit between the existing ones
#'
#' @param path Filepath to the panel .csv, or a data.frame with column
#'   Fluorophore
#' @param NumberDetectors Number of detectors for the instrument to pull
#'   references
#' @param TheCutoff Default is 0.9, determines HighOverlap count threshold
#' @param returnAll Whether to return all variants, default is FALSE
#' @param filename Desired name for the output .csv
#' @param outpath Path Location to store the the output .csv
#' @param returnCSV Whether to return as a csv to designated outpath and
#'   filename, a TRUE/FALSE
#'
#' @importFrom utils read.csv write.csv
#' @importFrom dplyr arrange bind_rows desc filter group_by pull relocate
#'   select slice ungroup
#' @importFrom purrr map
#'
#' @return A csv containing selected fluorophores
#' @export
#'
#' @examples
#' Folder_Location <- system.file("extdata", package = "Luciernaga")
#' ThePanelLocation <- list.files(Folder_Location, pattern = "^Panel.csv",
#'  full.names = TRUE)
#' OutPath <- file.path("C:", "Users", "JohnDoe", "Desktop")
#'
#' ProspectiveAdditions <- QC_ProspectiveAdditions(path = ThePanelLocation,
#' NumberDetectors = 64, TheCutoff = 0.9, returnAll = FALSE,
#' returnCSV = FALSE, filename = "ProspectiveAdditions", outpath = OutPath)
#'
QC_ProspectiveAdditions <- function(path, NumberDetectors, TheCutoff = 0.9,
                                     returnAll = FALSE, filename, outpath,
                                     returnCSV) {
  ReferenceData <- InstrumentReferences(NumberDetectors = NumberDetectors)

  if (!is.data.frame(path)) {
    TheList <- read.csv(path, check.names = FALSE)
    TheList <- TheList |> pull(Fluorophore)
  } else {
    TheList <- path |> pull(Fluorophore)
  }
  #x <- TheList[1]

  TheReferenceList <- ReferenceData |> dplyr::filter(Fluorophore %in% TheList)

  TheDetectorContained <- TheReferenceList |> group_by(Fluorophore) |>
    arrange(desc(AdjustedY)) |> slice(1) |> select(Fluorophore, Detector) |>
    ungroup()

  TheOccupiedDetectors <- TheDetectorContained |> pull(Detector) |> sort()

  TheOtherFluors <- ReferenceData |> dplyr::filter(!Fluorophore %in% TheList)
  TheOtherDetectors <- TheOtherFluors |> group_by(Fluorophore) |>
    arrange(desc(AdjustedY)) |> slice(1) |> select(Fluorophore, Detector) |>
    ungroup()
  TheOtherDetectors <- TheOtherDetectors |> arrange(Detector)
  TheOtherDetectors <- TheOtherDetectors |>
    dplyr::filter(!Detector %in% TheOccupiedDetectors)

  TheOtherDetectors <- TheOtherDetectors |>
    relocate(Detector, .before = Fluorophore)

  PossibleDetectors <- TheOtherDetectors |> pull(Detector) |>
    unique() #%>% length()

  # x <- PossibleDetectors[1]
  DataUnlocked <- map(.x = PossibleDetectors, .f = Comparison,
    TheOtherDetectors = TheOtherDetectors, TheList = TheList,
    ReferenceData = ReferenceData, TheCutoff = TheCutoff) |> bind_rows()

  DataUnlocked$RankValue <- as.numeric(DataUnlocked$RankValue)

  if (returnAll == TRUE) {
    PossibleLocations <- DataUnlocked |>
      arrange(TheDetector, RankValue)
  } else {
    PossibleLocations <- DataUnlocked |> group_by(TheDetector) |>
      arrange(RankValue) |> slice(1) |> ungroup()
    TheOutput <- PossibleLocations |> arrange(RankValue)
  }

  if (returnCSV == TRUE) {
    TheFileName <- paste0(filename, ".csv")
    StoreHere <- file.path(outpath, TheFileName)
    write.csv(TheOutput, StoreHere, row.names = FALSE)
  }

  return(TheOutput)
}